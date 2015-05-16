#!/bin/sh

python /home/vachasp2/HGT-Simulations/concatenate.py $1/all-genes.phylip > $1/concat.fasta

x=$1/concat.fasta

o=$1/concatenatedtree.genes1000
test -s $o && exit 1

$HOME/bin/fasttree -nt -gtr  -nopr  $x > $o

rm $x.gz
gzip $x

