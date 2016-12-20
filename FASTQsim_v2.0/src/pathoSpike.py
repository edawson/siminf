'''
This class spikes in silico reads from provided source files into a background FASTQ file.
 @author        Anna Shcherbina (mailto: anna.shcherbina@ll.mit.edu)
License:          GNU GPL license (http://www.gnu.org/licenses/gpl.html)  

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''

import sys 
import io 
import os 
from Params import *
from random import *
from helpers import *  
import numpy
from MultiprocessingHelpers import * 
from multiprocessing import Process, Queue, current_process,freeze_support,Manager

def main():    
    print "Parsing input arguments...\n"
    if len(sys.argv)<5:
        raise SyntaxError("usage: python pathoSpike.py [background fastq or -nobackground] [-params <parameters.csv> OR -platform <target_platform, one of \"ion\", \"illumina\",\"pacbio\",\"roche\">] [<coverage> <pathogen.fasta>] [-o output prefix] [-threads number of threads, default=1] [-plothistogram]")
    #determine whether to perform plotting 
    plothist=False
    plotCollection=None
    if sys.argv.__contains__('-plothistogram'): 
        plothist=True
        import matplotlib.pyplot as plt
        from plotSpiked import *
        plotCollection=plotSpiked();
        plotindex=sys.argv.index('-plothistogram') 
        sys.argv.pop(plotindex) 

    #determine whether we want to spike into a background file
    usebackground=True 
    if sys.argv.__contains__('-nobackground'): 
        usebackground=False
        backindex=sys.argv.index('-nobackground') 
        sys.argv.pop(backindex) 
    if usebackground:
        background_reads = sys.argv[firstarg]
        print "background_reads: " + background_reads + "\n"
    else:
        background_reads="NotSpecified"
        print "no background fastq file will be used for this run\n" 
    param_index=-1 
    platform_index=-1
    if "-params" in sys.argv: 
        param_index=sys.argv.index("-params") 
    if "-platform" in sys.argv: 
        platform_index=sys.argv.index("-platform") 
    if param_index==-1: 
        platform_index=sys.argv.index("-platform") 
        if platform_index==-1: 
            raise SyntaxError("No -params file or -platform designation specified. Please provide a paramater file or indicate the platform to use with default parameters\nusage: python pathoSpike.py [background fastq or -nobackground] [-params <parameters.csv> OR -platform <target_platform>] [<coverage> <pathogen.fasta>] [-o output prefix] [-threads number of threads, default=1] [-plothistogram]")
    if platform_index > -1:
        platform=sys.argv[platform_index+1] 
        if platform=='ion': 
            parameter_file_set="params/ion/ion_summary.csv" 
        elif platform=='illumina': 
            parameter_file_set="params/illumina/illumina_summary.csv" 
        elif platform=='roche': 
            parameter_file_set="params/roche/roche_summary.csv" 
        elif platform=='pacbio': 
            parameter_file_set="params/pacbio/pacbio_summary.csv" 
        else: 
            raise SyntaxError("Invalid platform name specified. Must be on of \"ion\" (IonTorrent), \"illumina\", \"pacbio\", \"roche\"")

    else:
        parameter_file_set = sys.argv[param_index+1]
    print "parameter_file_set: " + str(parameter_file_set) + "\n"
    

    #check to see if user specified an output file. 
    if "-o" in sys.argv: 
        outindex=sys.argv.index('-o') 
        outputfilename=sys.argv[outindex+1]
        sys.argv.pop(outindex+1)
        sys.argv.pop(outindex)

    else: 
        outputfilename= background_reads.split('.')[0]


    #check to see the number of threads to use 
    if "-threads" in sys.argv: 
        threadindex=sys.argv.index('-threads') 
        numthreads=int(sys.argv[threadindex+1]) 
        sys.argv.pop(threadindex+1) 
        sys.argv.pop(threadindex) 
    else: 
        numthreads=1 

    #may have multiple sources of pathogen DNA (i.e. example uses E.coli and plasmid) 
    pathogen_coverage_presplit = [] 
    pathogen_source_presplit = [] 
    is_circular_presplit = [] 
    
    for i in range(len(sys.argv)): 
        if sys.argv[i].__contains__('-source'): 
            pathogen_coverage_presplit.append(float(sys.argv[i+1])) 
            pathogen_source_presplit.append(sys.argv[i+2]) 
            is_circular_presplit.append(sys.argv[i+3]) 
    if len(pathogen_source_presplit)==0: 
        raise SyntaxError("No source fastq files were specified for performing the spiking operation. Please indicate at least one source FASTQ file using the format -source <coverage level> <FASTQ file> <circular genome?>") 
    print "paramter_file_set:"+str(parameter_file_set) 
    params = Params(8,parameter_file_set,background_reads,plotCollection)
    print "obtained params"
    readSource = []  
    revcompSource = [] #the reverse complement of the source files. 
    pathogenSourceLength = []
    is_circular=[] 
    pathogen_coverage=[] 
    pathogen_source=[] 
    for i in range(len(pathogen_source_presplit)):
        
        readSrcF = open(pathogen_source_presplit[i],'r').read().replace('\r\n','\n').replace('\r','\n')
        readSrcF=readSrcF.split('>')
        while '' in readSrcF:
            readSrcF.remove('')
        for line in readSrcF:
        #pull out the sequences from the FASTA input file and add them to the bases to be spiked in                                                                                                                                         
            firstnewline=line.find('\n')
            line_header=line[0:firstnewline]
            line=line[firstnewline+1::].replace('\n','')
            readSource.append(line)
            revcompSource.append(getRevComp(line))
            pathogenSourceLength.append(len(line))
            is_circular.append(is_circular_presplit[i])
            pathogen_coverage.append(pathogen_coverage_presplit[i])
            pathogen_source.append(pathogen_source_presplit[i]+"_"+line_header)
    spikedReads = []

#If -plothistogram is specified, keep track of actual position/size for insertions/deletions/mutations to ensure that the proper distribution/ parameters are being used to generate these values. This step should be ommitted if greater speed is desired. 
    if plothist:
        actualReadLengths=dict() 
        actualDelPos = dict() 
        actualInsertPos = dict() 
        actualMutPos = dict()         
        actualDelSize = dict() 
        actualInsertSize = dict()
        actualPosCount = dict() 

    task_queue=Queue() 
    done_queue=Queue() 

    print "starting worker processes"
    #start worker processes 
    procs=[] 
    for i in range(numthreads): 
        p=Process(target=read_gen_worker,args=(task_queue,done_queue,params))
        p.daemon=True
        p.start() 
        procs.append(p) 

    
    print "generating reads\n"
    proc_return_data =dict()
    for i in range(len(readSource)):
        proc_return_data[i]=0
        src = readSource[i]
        revComp = revcompSource[i]
        circularity = bool(is_circular[i]) 
        coverage = pathogen_coverage[i] 
        basesToGenerate = round(coverage*len(src)) 
        task_queue.put((i,src,revComp,circularity,coverage,basesToGenerate,plothist))
        
    
    #get items from the output queue 
    totalInputs=len(readSource) 
    returned = 0 
    while True: 
        if returned % 10000: 
            print "generated spike-in reads for " + str(returned) + " inputs, out of "+str(totalInputs) 
        if done_queue.empty()==False:
            [src_index,posCountDict_instance,spikedReads_instance]=done_queue.get() 
            if (posCountDict_instance==None) and (spikedReads_instance==None): #failure! 
                #try to reprocess this src up to 5 times 
                proc_return_data[src_index]+=1 
                if proc_return_data[src_index]>4: 
                    #failed to process this block 
                    returned+=1 
                    print "WARNING: 5 failed attempts to generate reads for input source " + str(src_index)+ "; aborting\n" 
                else: 
                    src = readSource[src_index]
                    revComp = revcompSource[src_index]
                    circularity = bool(is_circular[src_index]) 
                    coverage = pathogen_coverage[src_index] 
                    basesToGenerate = round(coverage*len(src)) 
                    task_queue.put((src_index,src,revComp,circularity,coverage,basesToGenerate,plothist))      
            else:
                if plothist: 
                    actualPosCount=mergeDict(actualPosCount,posCountDict_instance)
                spikedReads=spikedReads+spikedReads_instance
                returned+=1
        elif returned >=totalInputs: 
            break; 

    #add a 'STOP' to each task queue to indicate that no more tasks will be passed in 
    for i in range(numthreads): 
        task_queue.put('STOP')             
            
            
    allspiked=0
    for instance in spikedReads: 
        allspiked+=instance[1]
    print "BASES GENERATED:"+str(allspiked)+'\n' 

    
    del readSource # we no longer need this variable. 
    del revcompSource
    

    #make sure the daemon processes are dead 
    for p in procs: 
        p.terminate() 


    print "processing reads (generating insertions, mutations, deletions)...\n" 
    processedReads = [] 
    processedQuality=[] 
    processedSource=[] 
    processedStart=[] 
    processedEnd=[] 

    #dispatch each read to be processed in parallel. 
    task_queue=Queue() 
    done_queue=Queue() 
    #start worker processes 
    procs=[] 
    for i in range(numthreads): 
        p=Process(target=read_mutate_worker,args=(task_queue,done_queue,params))
        p.daemon=True
        p.start() 
        procs.append(p) 
    putval=0;
    proc_return_data=dict(); 
    for i in range(len(spikedReads)): 
        proc_return_data[i]=0
        entry = spikedReads[i] 
        task_queue.put((entry,i,plothist))
        putval+=1 
    returned = 0; 
    oldval=0; 
    while True:
        if returned %10000==0: 
            print "processsed:"+str(returned)+' spiked reads out of ' + str(len(spikedReads))+'\n' 
        if done_queue.empty()==False: 
            [src_index,processed_read,processed_qval,stat_dict,s,e,read_source]=done_queue.get() 
            if (processed_read == None) and (processed_qval==None) and (stat_dict==None):
                proc_return_data[src_index]+=1 
                if proc_return_data[src_index]>4: 
                    returned+=1 
                    print "WARNING: 5 failed attempts to mutate read:"+spiked_reads[src_index]+"; aborting\n"
                else:
                    entry=spikedReads[src_index] 
                    task_queue.put((entry,src_index,plothist))
            else:
                returned+=1;
                processedReads.append(processed_read)
                processedQuality.append(processed_qval) 
                processedSource.append(pathogen_source[read_source]) 
                processedStart.append(s) 
                processedEnd.append(e) 
                if plothist:
                    #update all relevant dictionaries 
                    if 'readLength' in stat_dict: 
                        if stat_dict['readLength'] not in actualReadLengths: 
                            actualReadLengths[stat_dict['readLength']]=1 
                        else: 
                            actualReadLengths[stat_dict['readLength']]+=1 
                    if 'posCount' in stat_dict: 
                        for posval in stat_dict['posCount']: 
                            if posval not in actualPosCount: 
                                actualPosCount[posval] = 1
                            else: 
                                actualPosCount[posval]+=1 

                    if 'mutPos' in stat_dict: 
                        if stat_dict['mutPos'] not in actualMutPos: 
                            actualMutPos[stat_dict['mutPos']]=1 
                        else:
                            actualMutPos[stat_dict['mutPos']]+=1 
                    if 'insertPos' in stat_dict: 
                        if stat_dict['insertPos'] not in actualInsertPos: 
                            actualInsertPos[stat_dict['insertPos']]=1 
                        else: 
                            actualInsertPos[stat_dict['insertPos']]+=1 
                    if 'insertSize' in stat_dict: 
                        if stat_dict['insertSize'] not in actualInsertSize: 
                            actualInsertSize[stat_dict['insertSize']]=[1,0]  
                        else: 
                            actualInsertSize[stat_dict['insertSize']][0]+=1 
                    if 'insertRepeat' in stat_dict: 
                        if stat_dict['insertRepeat'] not in actualInsertSize: 
                            actualInsertSize[stat_dict['insertRepeat']]=[0,1] 
                        else: 
                            actualInsertSize[stat_dict['insertRepeat']][1]+=1 
                    if 'delPos' in stat_dict: 
                        if stat_dict['delPos'] not in actualDelPos: 
                            actualDelPos[stat_dict['delPos']]=1 
                        else: 
                            actualDelPos[stat_dict['delPos']]+=1 
                    if 'delSize' in stat_dict: 
                        if stat_dict['delSize'] not in actualDelSize: 
                            actualDelSize[stat_dict['delSize']]=[1,0] 
                        else: 
                            actualDelSize[stat_dict['delSize']][0]+=1 
                    if 'delRepeat' in stat_dict: 
                        if stat_dict['delRepeat'] not in actualDelSize: 
                            actualDelSize[stat_dict['delRepeat']]=[0,1] 
                        else: 
                            actualDelSize[stat_dict['delRepeat']][1]+=1
        elif returned >=len(spikedReads): 
            print "finished mutating reads!"
            break; 
        else:
            sleep(1)
            if returned==oldval:                 
                print "processsed:"+str(returned)+' spiked reads out of ' + str(len(spikedReads))+'\n' 
            oldval=returned; 



   #add a 'STOP' to each task queue to indicate that no more tasks will be passed in 
    for i in range(numthreads): 
        task_queue.put('STOP') 
    #make sure the daemon processes are dead 
    for p in procs: 
        p.terminate() 

    del spikedReads 

    #plot statistics for readLength, insertions, deletions, and mutations 

    if plotCollection:
        plotCollection.actualReadLengths=actualReadLengths;         
        plotCollection.actualPosCount=actualPosCount; 
        plotCollection.actualDelPos=actualDelPos; 
        plotCollection.actualMutPos=actualMutPos; 
        plotCollection.actualInsertPos=actualInsertPos; 
        plotCollection.actualInsertSize=actualInsertSize; 
        plotCollection.actualDelSize=actualDelSize; 

    if usebackground: 
    #insert spiked reads into the background 
        print "reading in background read block\n"                   
    #read in the portion of the background reads to be spiked. 
        backgroundReadsf = open(background_reads,"r") 
        backgroundSize = os.stat(background_reads).st_size 
        maxToRead = 1048576/10 #.1 MB chunk 
        numIterations = math.ceil(float(backgroundSize)/maxToRead) 
        numToSpike = math.ceil(float(len(processedReads))/numIterations)
     
        allSpikedFastq = [] 
        allInsertedQual=[] 
        allInsertedFasta = [] 
        allInsertedNames=[] 
    
    #parallelization code for read spiking 
        task_queue=Queue() 
        done_queue=Queue() 
    
        procs=[] 
        for i in range(numthreads): 
            p=Process(target=insert_spiked_reads_worker,args=(task_queue,done_queue))
            p.daemon=True 
            p.start() 
            procs.append(p) 

        proc_return_data=dict(); 
        for iterval in range(int(numIterations)):
            backgroundReads= backgroundReadsf.read(maxToRead) 
            backgroundReads=backgroundReads.replace('\r','') 
        #make sure we end at a complete FASTQ sequence. 
            if '\n' not in backgroundReads: 
                readmany=True
            else:
                lastnewline=backgroundReads[::-1].index('\n') 
                if lastnewline==0:
                    readmany=True
                elif backgroundReads[::-1][lastnewline-1:lastnewline+1]=='+\n': 
                    readmany=False
                else: 
                    readmany=True 
            if readmany: 
                while backgroundReads.endswith('+\n')==False: 
                    newchar=backgroundReadsf.read(1)
                    if newchar=="":
                        break; 
                    if newchar!='\r':
                        backgroundReads+=newchar 
            while backgroundReads.endswith('\n')==False: 
                newchar=backgroundReadsf.read(1) 
                if newchar=="": 
                    break; 
                if newchar!='\r': 
                    backgroundReads+=newchar 
            backgroundReads = backgroundReads.split("\n") 

            #get the fastq reads from this file 
            qualwithname,qualDict = fastqextract(backgroundReads)

            #pick a subset of the reads ( randomly) to spike into this portion of the file. 
            subset_processed=dict()
            longest_allowed = max(qualDict.keys()) 

            for i in range(int(numToSpike)): 
                if len(processedReads) > 0: 
                    pos = random.randint(0,len(processedReads)-1) 
                    new_read = processedReads.pop(pos)
                    new_qual=processedQuality.pop(pos); #get the matching quality score. 
                    new_source_id=processedSource.pop(pos); 
                    new_start_index=processedStart.pop(pos); 
                    new_end_index=processedEnd.pop(pos); 
                    
                    if len(new_read)> longest_allowed: 
                        new_read=new_read[0:longest_allowed] 
                        new_qual=new_qual[0:longest_allowed]
                    if subset_processed.__contains__(len(new_read)): 
                        subset_processed[len(new_read)].append(([new_read,new_qual,new_source_id,new_start_index,new_end_index])) 
                    else: 
                        subset_processed[len(new_read)]=[([new_read,new_qual,new_source_id,new_start_index,new_end_index])]  

            #print "inserting reads into background slice " + str(iterval)+"\n" 
            #print "starting insertSpikedreads\n"
            proc_return_data[iterval]=[0,backgroundReads,subset_processed,qualwithname,qualDict]
            task_queue.put((iterval,backgroundReads,subset_processed,qualwithname,qualDict,pathogen_source))


        returned =0;
        while True: 
            print "spiked "+str(returned) + " file blocks, out of "+str(numIterations)+"\n" 
            if done_queue.empty()==False: 
                [src_index,spikedFastq,insertedQual,insertedFasta,insertedReadNames]= done_queue.get() 
                if (spikedFastq==None) and (insertedQual==None) and (insertedFasta==None) and (insertedReadNames==None):
                    proc_return_data[src_index][0]+=1 
                    if proc_return_data[src_index][0]>4: 
                        returned+=1 
                        print "WARNING: 5 failed attempts to spike background read block # "+str(src_index)+"; aborting\n" 
                    else: 
                        task_queue.put((src_index,proc_return_data[src_index][1],proc_return_data[src_index][2],proc_return_data[src_index][3],proc_return_data[src_index][4]))
                else:
                    allSpikedFastq=allSpikedFastq + spikedFastq
                    allInsertedQual = allInsertedQual+insertedQual 
                    allInsertedFasta = allInsertedFasta+insertedFasta 
                    allInsertedNames = allInsertedNames+insertedReadNames
                    returned+=1; 
                    print "spiked "+str(returned) + " file blocks, out of " + str(numIterations)+'\n' 
            elif (task_queue.empty()==False) or (returned < int(numIterations)):
                sleep(1) 
            else: 
                break;     
        backgroundReadsf.close()

        #add a 'STOP' to each task queue to indicate that no more tasks will be passed in 
        for i in range(numthreads): 
            task_queue.put('STOP') 

        for p in procs: 
            p.terminate() 
        print "writing spiked reads to separate file\n"
    else:
        allInsertedNames=[]# ['@'+outputfilename+"_"+str(i) for i in range(len(processedReads))]
        for var in range(len(processedSource)): 
            name="@"+processedSource[var].split('/')[-1]+"|"+str(processedStart[var])+"|"+str(processedEnd[var])
            allInsertedNames.append(name) 
        allInsertedFasta=processedReads; 
        allInsertedQual=processedQuality; 

    print "writing spiked output file"
    spikedOutputFastaf = open(outputfilename+"_Spiked.fq","w")     
    spikedreadlengthdict=dict() 
    for pindex in range(len(allInsertedFasta)):
        p = allInsertedFasta[pindex] 
        if spikedreadlengthdict.__contains__(len(p)):
            spikedreadlengthdict[len(p)]+=1 
        else: 
            spikedreadlengthdict[len(p)]=1 
        q = allInsertedQual[pindex]
        name = allInsertedNames[pindex] 
        spikedOutputFastaf.write(name+"\n")
        spikedOutputFastaf.write(p.upper()+"\n") 
        spikedOutputFastaf.write("+\n") 
        spikedOutputFastaf.write(q.upper()+'\n')
    spikedOutputFastaf.close() 

    
    if usebackground:
        fullOutputFastaf = open(outputfilename+"_Full.fq","w")
        for f in allSpikedFastq: 
            fullOutputFastaf.write(f.upper()+"\n") 
        fullOutputFastaf.close()    

    #generate summary graphs for background distributions, regression curves, and spiked reads 
    if plotCollection:
        plotCollection.plotAll(); 
if __name__ == '__main__':
    main() 
