mkdir ref 
cd ref
mkdir human 
cd human 
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT 
do                                                                                                                                                                                                                                          
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38_chr$i.fa.gz                                                                                                                                       gunzip *.fa.gz
done
for i in alts unlocalized unplaced                                                                                                                                                                                                       
do                                                                                                                                                                                                                                          
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38_$i.fa.gz                                                                                                                                          
gunzip *.fa.gz
done     
cd ..                                                                                                                                                                   
cat human/*.fa > hs_ref_GRCh38.fa
cd ..


