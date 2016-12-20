cd burke_seqs/
for i in `ls`
do
    cd $i
    jname=`basename $(ls | grep ".fsa$" | head -n 1 | cut -f 3 -d "_" ) .fsa`
    echo "Seq name: $jname"
    cp $(ls | grep ".fsa$" | head -n 1) example_${jname}.fasta
    sed -i "s/>.*/>$jname/g" example_${jname}.fasta
    ~/Dropbox/siminf/bed_extract.sh $(ls | grep ".fsa$" | head -n 1) && mv tmp example_${jname}.bed
    cd ..
done
cd ../
