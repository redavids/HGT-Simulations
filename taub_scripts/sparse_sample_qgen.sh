#!/bin/bash

quartetcounts=${4}

if [ ! -f $quartetcounts ]
then
    quartets/quartet-controller.sh $2 $quartetcounts
fi

python quartetsparsesample-driver.py $1 $quartetcounts $3
