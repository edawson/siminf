## Generating HPV coinfections and detecting them with a linear learner

## Generate coinfections
Create a mix.txt file that contains your coinfections / pure samples.
It will look like this:
        
            # Comments are fine in these files
            # here's a pure sample:
            A4  1.0 3000    /path/to/simreads.fq

            # here's a coinfection
            A1  0.5 3000    /path/to/a1_reads.fq
            A4  0.5 3000    /path/to/a4_reads.fq

Next, we'll make these files with the make\_mix.py script.
Our `-i` parameter is our mix.txt and our `-l` is the approximate genome length.

        python make_mix.py -i mix.txt -l 7000

This will call seqtk the correct way a bunch of times, concatenate all of our sub-infection
files together, and remove any intermediate files. You'll get an output file with a name like "D3\_0-0.491197832028\_C\_0-0.0175280259015\_B1\_0-0.014299287417\_D1\_0-0.0127703817737\_D2\_0-0.464204472879\_7898\_.fq"

with the pattern [label\_proportion]\*\_RandInt\_.fq.

## Train your model

To train your model, first classify each read using RKMH. We'll call our file containing reads "readsfile.fq", but
our real input would be the file with the long name that came out of our simulator above. For the references, use
the reference genomes from which the reads originally came, all in one fasta file. We'll use a kmer size of 20,
a sketch size of 4000, and we'll flag any reads which have fewer than 20 kmers between their best and second-best classifications:
            
            rkmh classify -f readsfile.fq -r refs.fa -k 20 -s 4000 -D 20 > rkmh.cls.txt

Now have to transform our initial result to VW format. we can do this with a script in the scripts dir. We'll also
label our infection status here using the `-C` flag - 0 for not coinfected, 1 for coinfected:
        
        python scripts/vwize.py -i rkmh.cls.txt -C 1 >> vecs.txt

This will create **just one line per input file**, essentially compressing our initial data waaaaay down. We'll need to run this across many, many files
to generate a sufficient training set (ideally on the order of 10,000).

We could train off this file, but it might be helpful to collapse subtypes (lineages starting with the same letter) to a single feature:

        cat vecs.txt | python collapse_subtypes.py > collapsed.vecs.txt

Now we can use this file to train our learner. We'll use vowpal-wabbit, mostly because it's easy and plenty good for most purposes.
Usually we would generate a holdout set but in this case we'll train on a random subset of the whole:

        cat collapsed.vecs.txt | shuf | head -n 1000 | vw --passes=20 --cache_file cache.f --binary --interactions vvv -f trained.model

This will generate a trained model called "trained.model". We allow vw to recycle data in an attempt to improve learning (--passes=20). We also
tell it to use a binary error model (--binary), use a cache file for the passes called "cache.f", and to allow cubic interactions between features.
Interactions can be though of as combinations of features; if we had two features [A,B], interacting them with quadratic interactions (i.e. 2 layers) would
allow the model to generate a feature (A\*B), extending our model to have features [A, B, (A\*B)]. This is very useful for coinfections since we expect
that our features (strain proportions) should influence one another.

## Apply (test) your model

To apply our model, we will use our full dataset and print results to stdout:

        cat collapsed.vecs.txt | shuf | vw -i trained.model -p /dev/stdout | tee classified.txt

classified.txt will contain one line per coinfection/pure infection, with two columns. The first is the
logistic score from the model and the second one is the label we gave to vw. We can quantify our performance using the conf\_mat.py script:

        python conf_mat.py -i classified.txt

So far, we're seeing low false positive rates but high false negative rates (lots of coinfections slip through). That's not good, but maybe our model will improve with
more data.
