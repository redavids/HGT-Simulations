#!/bin/sh

python /scratch/01721/smirarab/astral2/secondround/concatenate.py $1/all-genes.phylip $1/wellresolvedgenes $2 > $1/concat.fasta.$2

x=$1/concat.fasta.$2

o=$1/concatenatedtree.genes$2

test -s $o && exit 1

$HOME/bin/fasttree -nt -gtr -quiet -nopr  $x > $o

rm $x.gz
gzip $x

