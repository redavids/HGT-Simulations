#!/bin/bash

x=$1/all-genes.phylip

gunzip $x.gz

o=${x/all-genes.phylip/estimatedgenetre}

$HOME/bin/fasttree -nt -gtr -quiet -nopr -gamma -n 1000 $x > $o

#-quiet means suppresses reporting info
# -nopr suppress progress indicator
# -gamma rescales tree length after CAT approximation 

gzip $x

wc -l $o
