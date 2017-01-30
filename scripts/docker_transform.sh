infile=$1
collapse=$2
if [ "$2" == "collapse" ]
    then /bin/python vwize.py -i $infile -C 1 | /bin/python collapse_subtypes.py
else
    python /bin/vwize.py -n -i $infile -C 1
fi

