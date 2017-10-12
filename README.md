siminf
------
For simulating HPV coinfections  
Eric T Dawson  
December 2016

### Requirements
To regenerate the data in this repo, you'll need to have BEDTools installed. You'll also need wgsim,
ART, or some other FASTQ simulator. I used ART because it's fast and good enough.


## Generating test data:
1. Simulate fastqs for each strain    

``./scripts/make_reads_single_strain.sh <path/to/strain/fastq.fq>``  

You might want to simulate a single fastq for each sublineage (i.e. combine the strain-level fastqs within a sublineage to build a composite reference fastq).
You can use the make_reads_single_lineage.sh script to do this; pass the directory containing a set of fastqs as its sole argument.

2. Make a two-column (tab separated) file that describes the strain label -> path/to/file relationships.
Your file should look like:  

            A1	/home/user/a1/bigA1.fq  
            B1	/home/user/b1/bigb1.fq  
            C1	/home/user/c1/bigc1.fq
  
We'll save it as `path_file.txt`. 

3. This file will become the input to a script that makes descriptions for simulated infections. To simulate 100 single-strain infections:  

`` python scripts/make_random_mix.py -n 100 -i path_file.txt > desc_file.txt``  
**NB:** We can also pass `-m <coverage>` to set a minimum required coverage for each file, and `-c ` to simulate coinfections rather than single-strain ones.

4. We'll use our new `desc_file.txt` as an argument to sample reads from our fastqs in the right proportions to generate our desired coinfections.
We also pass an approximate genome length to help calculate the number of reads needed to achieve our desired coverage.
**Warning: This step takes a long time, hammers your drive with lots of write and can make a lot of data very quickly**  

``python scripts/make_mix.py -i desc_file.txt -l 7000 ``


## Data description
### 1: HPV subtype coinfection, simulation, Ion Torrent
   Instructions:
        1. Generate large sets of simulated reads for each subtype
        using ART (MiSeq v3 profile) and FASTQSim (Ion Torrent profile)
        2. Subsample the simulated fastqs for the desired strains using seqtk.
        Make sure that the total number of reads matches the desired coverage
        and that the number of each from each subtype matches the desired proportion.

    Available data:
        [ ] 0.15_A4-0.15_D3-0.15_A1-0.15_B1-0.4_C.fq : 15% each of A4, A1, B1, D3, and 40% C, 30X coverage
        [ ] 0.25_A4-.75_A1.fq : 25% A4 and 75% A1, 30X coverage
        [ ] 0.01_A1-0.01_D1-0.98_A4 : 1% A1, 1% B1, 98% D1, 10X, 30X and 100X
        [ ] 0.01_A1-0.01_A4-0.20_B-0.78_D1 : 1/% A1, 1% A4, 20% B and 78% D1 at 10X, 30X, and 100X
### 2: HPV subtype coinfection, simulation, Nanopore
    Instructions:
        1. Simualte large sets of simulated reads for each subtype,
        this time using wgsim with high error/indel rates (25% / 8%) and long reads (6.5kb average),
        which is approximately what we'd expect from bad 1D nanopore reads before correction.

    Data Availability:
        [ ] 0.15_A4-0.15_D3-0.15_A1-0.15_B1-0.4_C.fq : 15% each of A4, A1, B1, D3, and 40% C, 1000X coverage
        [ ] 0.25_A4-.75_A1.fq : 25% A4 and 75% A1, 1000X coverage
        [ ] 0.01_A1-0.01_D1-0.98_A4 : 1% A1, 1% D1, 98% D1, 100X, 1000X and 5000X
        [ ] 
### 3: HPV subtype coinfection, simulation, Ion Torrent, amplicon restricted
    Instructions:
        1. For each subtype, restrict its reference sequence to amplicon regions.
        2. Generate reads as in the initial Ion Torrent simulation, using ART and FASTQSim.
        3. Mix reads by subsampling and combining reads at the desired ratios using seqtk.

    Data Availability:
        [ ]
        [ ]
        [ ]
        [ ]
### 4: HPV subtype coinfection, real, Ion Torrent
    Instructions:
        1. Based on manual annotation and previous pipeline runs, locate samples which
        appear to be homogeneous.
        2. Subsample these BAMs/FASTQs using seqtk
        3. Combine at the desired ratios and for the desired coverage.

    Data Availability:
        [ ]
        [ ]
        [ ]
        [ ]
### 5: HPV subtype coinfection, multiple HPV types
    Instructions:
        1. Generate fastqs for HPV types (in addition to previously generated subtypes).
        2. Subsample FASTQs using seqtk, mixing both types and subtypes in the desired
        proportions and final coverage.

    Data Availability:
        [ ]
        [ ]
        [ ]
        [ ]
