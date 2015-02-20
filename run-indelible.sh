#!/bin/bash

cd $1

/projects/sate7/tools/bin/indelible

cat *phy | sed '/^ *$/d' > all-genes.phylip

rm *phy
rm *fas
