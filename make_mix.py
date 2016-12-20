import argparse
import subprocess
import math
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


## Returns an temp file name to be used later.
def make_sample(sample, percent_comp, num_reads, samp_file):
    ## First: generate temp file name
    t_file = "_".join([sample, str(percent_comp), ".tmp"])
    cmd = "seqtk/seqtk sample -s42 " + samp_file + " " + str(num_reads) + " > " + t_file
    print cmd
    subprocess.call(cmd, shell=True)
    return t_file

if __name__ == "__main__":
    args = parse_args();

    mid_files = []

    cov = -1
    with open(args.sample_file, "r") as ifi:
        for line in ifi:
            if line.startswith("#"):
                continue
            tokens = line.strip().split("\t")
            ## use empty lines as file delims
            ## collapse our N intermeidate files
            ## into an output file, clear the list of current intermed files,
            ## and write our output.

            if (len(tokens) == 1):
                outname = "_".join( [ "-".join([x.split(".")[0].strip("_"), "0." + x.split(".")[1].strip("_").strip("_")]) for x in mid_files ] ) + "_" + str(cov) + "_" + ".fq"
                for i in mid_files:
                    print "Writing i to: ", outname
                    ## just cat the mid_files to output
                    o_cmd = "cat " + i + " >> " + outname
                    subprocess.call(o_cmd, shell=True)
                    ## remove the mid_file
                mid_files = []
                cov = -1
            else:
                mid_files.append(make_sample(tokens[0], tokens[1], get_num_reads(tokens[2], args.genome_length, tokens[1]), tokens[3]))
                cov = int(tokens[2])
            

