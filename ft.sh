#!/bin/bash

x=$1/all-genes.phylip

gunzip $x.gz

o=${x/all-genes.phylip/estimatedgenetre}

$HOME/bin/fasttree -nt -gtr -quiet -nopr -gamma -n 1000 $x > $o

gzip $x

wc -l $o
