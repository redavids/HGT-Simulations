#!/bin/bash



quartetcounts=${3}

if [ ! -f $quartetcounts ]
then
    quartets/quartet-controller.sh $1 $quartetcounts
fi

cp $quartetcounts $2
