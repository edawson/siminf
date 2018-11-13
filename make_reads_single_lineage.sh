touch big.fa
rm big.fa
for i in `ls | grep fsa`
do
    cat $i >> big.fa
done
~/sandbox/siminf/art_bin_MountRainier/art_illumina -ss MSv3 -sam -i big.fa -l 150 -o `basename big.fa .fa`_$1 -l 250 -c 280000 -rs 42
