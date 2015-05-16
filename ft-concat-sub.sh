#!/bin/sh

python /home/vachasp2/HGT-Simulations/concatenate.py $1/all-genes.phylipm$3 <(seq 1 1000) $2 > $1/concat.fastam$3.$2

python /home/vachasp2/HGT-Simulations/concatenate.py $1/all-genes.phylipm$3 <(seq 1 1000) $2 > $1/concat.fastam$3.$2

x=$1/concat.fastam$3.$2

#o=$1/concatenatedtreem.genes$2
o=$4

test -s $o && exit 1

$HOME/bin/fasttree -nt -gtr  -nopr  $x > $o

rm $x.gz
gzip $x

