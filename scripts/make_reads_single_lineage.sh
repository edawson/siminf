touch big.fa
rm big.fa
for i in `find $1 | grep fsa`
do
    cat $i >> big_$(basename $1).fa
done
art_illumina -ss MSv3 -na -i big_$(basename $1).fa -o big_$(basename $1) -l 250 -f 50000 -rs 42
