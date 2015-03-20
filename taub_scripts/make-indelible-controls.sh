#!/bin/bash

for w in 50 100 150; do 
    for n in 0  500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000; do 
        for i in {1..50}; do
            #cd /Users/ruthdavidson/code/hgt-data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n/$i/
            #python ~/code/HGT-Simulations/taub_scripts/splitter.py shuftruegenetrees 
            python ~/code/HGT-Simulations/taub_scripts/splitter.py /Users/ruthdavidson/code/hgt-data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n/$i/ shuftruegenetrees 
            #cd /home/redavid2/HGT-Simulations/taub_scripts/data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n/$i/
            #unset j
            #for x in $(<shuftruegenetrees); do echo "$x" > shuftruegenetrees.$((++j)).trees; done
            #perl /Users/ruthdavidson/code/post_stidsim.pl /Users/ruthdavidson/code/hgt-data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n/$i 1
            #perl /Users/ruthdavidson/code/HGT-Simulations/post_stidsim.pl `pwd` 1
            #rm /Users/ruthdavidson/code/hgt-data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n/$i/g_trees*.trees
        done
    done
done

for w in 50 100 150; do 
    for n in 0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000; do 
        #perl /Users/ruthdavidson/code/post_stidsim.pl /Users/ruthdavidson/code/hgt-data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n 1
        perl /Users/ruthdavidson/code/HGT-Simulations/ruth_post_stidsim.pl /Users/ruthdavidson/code/hgt-data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n 1
    done
done
for w in 50 100 150; do 
    for n in 0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000; do 
        for i in {1..50}; do
            rm /Users/ruthdavidson/code/hgt-data/RuthPerGeneNormSimulatedDatasets.single.highway/hw$w/noise$n/$i/g_trees*.trees
        done
    done
done
