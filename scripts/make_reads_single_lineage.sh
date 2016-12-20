
touch big.fa
rm big.fa
for i in `ls $1 | grep fsa`
do
    cat $i >> big.fa
done
art_illumina -ss MSv3 -na -i big.fa -o `basename big.fa .fa`_$2 -l 250 -f 50000 -rs 42
