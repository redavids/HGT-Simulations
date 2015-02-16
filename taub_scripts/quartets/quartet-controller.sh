#!/bin/sh

quartetcounts=${2}

if [ ! -s $quartetcounts ]
then
    sh ~/Dropbox/Programs/phylogenetics/quartets/quartet_count.sh $1 | perl ~/Dropbox/Programs/phylogenetics/quartets/summarize_quartets_stdin.pl > $2
fi


#sh quartet_count.sh $1 | perl summarize_quartets_stdin.pl > $2
