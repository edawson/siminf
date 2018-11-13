#x_name=$(cat $1 | grep -o "^>[A-Za-z0-9|]*" | cut -f 3 -d "|")
x_name=$(cat $1 | grep -o "^>[A-Za-z0-9|]*" | grep -o "[A-Za-z0-9|]*" )
echo $x_name
cat ~/Dropbox/siminf/HPV16Ref_Mapped_hg19.unmerged-detail_renamed.bed | sed "s/HPV16Ref/${x_name}/g" > tmp.$(echo $x_name | sed -r 's/\|/_/g')
