import argparse
import random
import multiprocessing as mp
import subprocess
import math
import os
## Usage: make_mix.py -s sample_file -o outfile
## sample file format: StrainName   PercentageComposition   TotalCoverage    File

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--sample-file", type=str, required=True, dest="sample_file")
    parser.add_argument("-l", "--genome-length", type=float, required=True, dest="genome_length")
    return parser.parse_args()

def get_num_reads(total_cov, genome_len, percent_comp):
    n_tot = float(total_cov) * (float(genome_len) / float(250))
    n_reads = math.ceil( float(n_tot) * float(percent_comp) )
    return int(n_reads)

def func(task):
    # subprocess.Popen(task, stdout=subprocess.PIPE, shell=True, preexec_fn=os.setsid)
    try:
        subprocess.call(task, shell=True)
    except KeyboardInterrupt, Exception:
        raise KeyboardInterrupt
    return


## Returns an temp file name to be used later.
def make_sample(sample, percent_comp, num_reads, samp_file):
    ## First: generate temp file name
    t_file = "_".join([sample, str(percent_comp), str(num_reads), str(random.randint(1,1000000)), ".tmp"])
    cmd = "seqtk/seqtk sample -s42 " + samp_file + " " + str(num_reads) + " > " + t_file
    print cmd
    return t_file, cmd

if __name__ == "__main__":
    args = parse_args();

    mid_files = []

    pool = mp.Pool(4)

    cov = -1
    o_names = []
    with open(args.sample_file, "r") as ifi:
        for line in ifi:
            if line.startswith("#") or line.startswith("\n"):
                continue
            tokens = line.strip().split("\t")
            ## use empty lines as file delims
            ## collapse our N intermeidate files
            ## into an output file, clear the list of current intermed files,
            ## and write our output.

            mid_files.append(make_sample(tokens[0], tokens[1], get_num_reads(tokens[2], args.genome_length, tokens[1]), tokens[3]))
            cov = int(tokens[2])
            outname = tokens[0] + "_" + str(cov) + "_" + ".fq"
            o_names.append(outname)

    work = [ o[1] for o in mid_files ]
                
    tmps = [ x[0] for x in mid_files ]


    try:
        pool.map_async(func, work).get(999999999999999)
    except KeyboardInterrupt:
        pool.terminate()
        exit()
                ## remove the mid_file

