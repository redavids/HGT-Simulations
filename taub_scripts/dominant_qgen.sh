#!/bin/bash



quartetcounts=${4}

if [ ! -f $quartetcounts]
then
    quartets/quartet-controller.sh $2 $quartetcounts
fi

python quartetscores-driver-dominant.py $quartetcounts $2
