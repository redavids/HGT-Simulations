#!/bin/bash



quartetcounts=${3}

if [ ! -f $quartetcounts ]
then
    quartets/quartet-controller.sh $1 $quartetcounts
fi

python quartetscores-driver-dominant.py $quartetcounts $2
