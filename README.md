siminf
------
For simulating HPV coinfections  
Eric T Dawson  
December 2016

### Requirements
To regenerate the data in this repo, you'll need to have BEDTools installed. You'll also need ART, a modern python install, and seqtk built.

### Process
1. simulate reads using art from a FASTQ:
    ``./make_reads_single_fa.sh ref.fa my_output_name``

2. make a file describing the desired lineages, their ratio, their total coverage,
   etc. A given line should like this, and multiple consecutive lines are included in the same
   file until an empty line is hit:
   ``A1 0.25    3000    /path/to/file``

   The fields are: lineageLabel, proportion, total coverage, path to file.

3. Subsample reads to desired proportions with the desired python wrapper script.
    ``python make_mix.py -i mix.txt -l 7000 ``  
    where the `-i` parameter is the mixfile from step 2 and the `-l` parameter is the genome length.

## Data description
### 1: HPV subtype coinfection, simulation, Ion Torrent
   Overview:
        1. Generate large sets of simulated reads for each subtype
        using ART (MiSeq v3 profile) and FASTQSim (Ion Torrent profile)
        ```` ````
        2. Subsample the simulated fastqs for the desired strains using seqtk.
        Make sure that the total number of reads matches the desired coverage
        and that the number of each from each subtype matches the desired proportion.

    Available data:

        [ ] 0.15_A4-0.15_D3-0.15_A1-0.15_B1-0.4_C.fq : 15% each of A4, A1, B1, D3, and 40% C, 3000X coverage
        [ ] 0.25_A4-.75_A1.fq : 25% A4 and 75% A1, 3000X coverage
        [ ] 0.01_A1-0.01_D1-0.98_A4 : 1% A1, 1% B1, 98% D1, 100X, 3000X and 10000X
        [ ] 0.01_A1-0.01_A4-0.20_B-0.78_D1 : 1/% A1, 1% A4, 20% B and 78% D1 at 100X, 3000X, and 10000X

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
