~/sandbox/siminf/art_bin_MountRainier/art_illumina -ss MSv3 -sam -i $1 -o `basename $1 .fa` -l 250 -f 50000
if [ ! -d alns ]; then mkdir alns; fi
if [ ! -d fq ]; then mkdir fq; fi
mv *.fq fq/
mv *aln alns
mv *sam alns
