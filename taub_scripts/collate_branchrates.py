from collections import defaultdict, namedtuple
import numpy as np
import os, sys

Run = namedtuple('Run', ['experiment', 'hgt', 'index', 'ngenes', 'filename', 'value'])

def parse_filename(folder, fname):
    parts = fname.split('_')
    value = open(folder + "/" + fname).read()
    if value:
        value = float(value)
    else:
        value = -1
#    parts[3] = int(parts[3])
#    if parts[3] == '.':
#        parts[3] = 1000
#    else:
#        parts[3] = int(parts[3].split('.')[1])
    return Run(parts[0], parts[1], parts[2], parts[3], fname, value) 
    
def sort_files(folder):
    runs = [parse_filename(folder, i) for i in os.listdir(folder)]
    runs = [i for i in runs if i.value >= 0]
    runs_averaged = []
    d = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
    for i in runs:
        d[i.experiment][i.hgt][i.ngenes].append(i.value)
    for experiment in sorted(d.keys()):
        print experiment
        d1 = d[experiment]
        for hgt in sorted(d1.keys()):
            d2 = d1[hgt]
            for ngenes in sorted([int(i) for i in d2.keys()]):
                runs_averaged.append( Run(experiment, hgt, 0, ngenes, '', np.average(d2[str(ngenes)])) )
    return d, runs_averaged

def printtable(methods_subset, hgtrates, genecounts, runs):
    for hgtrate in hgtrates:
        print "HGT ngenes", ' '.join([i.replace("missingbranchrate","") for i in methods_subset])
        for ngenes in sorted([int(i) for i in genecounts]):
            try:
                print str(hgtrate) + " " + str(ngenes) + " " + " ".join([['{:.1f}'.format(i.value * 100) for i in runs if i.hgt == hgtrate and i.ngenes == ngenes and i.experiment == m][0] for m in methods_subset]) 
            except IndexError:
                pass
                
        print
        print

def analyze(folder):
    d, runs = sort_files(folder)
    for experiment in sorted(d.keys()):
        print experiment
        d1 = d[experiment]
        for hgt in sorted(d1.keys()):
            d2 = d1[hgt]
            for ngenes in sorted([int(i) for i in d2.keys()]):
                print "|" + hgt + "\t|" + str(ngenes) + "\t|" + str(np.average(d2[str(ngenes)])) + "\t|" + str(np.std(d2[str(ngenes)])) + "\t|" + str(len(d2[str(ngenes)])) + "\t|"
        print 
    hgtrates = set([i.hgt for i in runs])
    methods = sorted(list(set([i.experiment for i in runs])))
    methods_true = sorted(list(set([i for i in methods if "true" in i])))
    methods_estimated = sorted(list(set([i for i in methods if "true" not in i])))
    methods_inv_true = sorted(list(set([i for i in methods_true if "invariants" in i])))
    methods_inv_true = [methods_inv_true[-1]] + methods_inv_true[:-1]
    methods_inv_estimated = sorted(list(set([i for i in methods_estimated if "invariants" in i])))

    methods_sparsesample_estimated = sorted(list(set([i for i in methods_estimated if "sparsesample" in i and  i[-1] =='p']))) + [i for i in methods if 'wqmc-estimated' in i]
    methods_sparsesample_true = sorted(list(set([i for i in methods_true if "sparsesample" in i and i[-1] =='p' ]))) + [i for i in methods if 'wqmc-true' in i]

    genecounts = set([i.ngenes for i in runs])


    printtable(methods_estimated, hgtrates, genecounts, runs)
    printtable(methods_true, hgtrates, genecounts, runs)

    printtable(methods_inv_estimated, hgtrates, genecounts, runs)
    printtable(methods_inv_true, hgtrates, genecounts, runs)

    printtable(methods_sparsesample_estimated, hgtrates, genecounts, runs)
    printtable(methods_sparsesample_true, hgtrates, genecounts, runs)

            

if __name__ == "__main__":
    folder = sys.argv[1]
    analyze(folder)
