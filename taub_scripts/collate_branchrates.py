from collections import defaultdict, namedtuple
import numpy as np
import os, sys

Run = namedtuple('Run', ['experiment', 'hgt', 'index', 'ngenes', 'filename', 'value'])

def add_astral(d):
    d['astral']['0'][   '10' ]= 0.130004 
    d['astral']['0'][   '25' ]= 0.086246 
    d['astral']['0'][   '50' ]=  0.05958 
    d['astral']['0'][  '100' ]= 0.044584 
    d['astral']['0'][  '200' ]= 0.029994 
    d['astral']['0'][  '400' ]= 0.024156 
    d['astral']['0'][ '1000' ]= 0.018326 
    d['astral']['02'][   '10' ]= 0.133754 
    d['astral']['02'][   '25' ]=    0.085 
    d['astral']['02'][   '50' ]= 0.057498 
    d['astral']['02'][  '100' ]= 0.034164 
    d['astral']['02'][  '200' ]= 0.027082 
    d['astral']['02'][  '400' ]= 0.020832 
    d['astral']['02'][ '1000' ]=  0.01666 
    d['astral']['05'][   '10' ]= 0.155414 
    d['astral']['05'][   '25' ]= 0.090002 
    d['astral']['05'][   '50' ]= 0.067084 
    d['astral']['05'][  '100' ]= 0.049162 
    d['astral']['05'][  '200' ]= 0.038748 
    d['astral']['05'][  '400' ]= 0.031246 
    d['astral']['05'][ '1000' ]=  0.02666 
    d['astral']['2'][   '10' ]= 0.148748 
    d['astral']['2'][   '25' ]= 0.092916 
    d['astral']['2'][   '50' ]= 0.072504 
    d['astral']['2'][  '100' ]=   0.0575 
    d['astral']['2'][  '200' ]= 0.045418 
    d['astral']['2'][  '400' ]=  0.03499 
    d['astral']['2'][ '1000' ]= 0.029572 

    d['astral']['20'][   '10' ]= 0.184164
    d['astral']['20'][   '25' ]=    0.106668
    d['astral']['20'][   '50' ]= 0.0675
    d['astral']['20'][  '100' ]= 0.05167
    d['astral']['20'][  '200' ]= 0.035414
    d['astral']['20'][  '400' ]= 0.03166
    d['astral']['20'][ '1000' ]=  0.023328
    d['astral']['50'][   '10' ]= 0.281252
    d['astral']['50'][   '25' ]= 0.157498
    d['astral']['50'][   '50' ]= 0.103336
    d['astral']['50'][  '100' ]= 0.072916
    d['astral']['50'][  '200' ]= 0.052496
    d['astral']['50'][  '400' ]= 0.034576
    d['astral']['50'][ '1000' ]=  0.020828

    d['astral-true']['0'][   '10' ]=  0.08375 
    d['astral-true']['0'][   '25' ]= 0.050412 
    d['astral-true']['0'][   '50' ]=  0.03416 
    d['astral-true']['0'][  '100' ]= 0.024998 
    d['astral-true']['0'][  '200' ]= 0.019988 
    d['astral-true']['0'][  '400' ]= 0.014154 
    d['astral-true']['0'][ '1000' ]= 0.007906 
    d['astral-true']['02'][   '10' ]= 0.100836 
    d['astral-true']['02'][   '25' ]= 0.059166 
    d['astral-true']['02'][   '50' ]= 0.039166 
    d['astral-true']['02'][  '100' ]= 0.027488 
    d['astral-true']['02'][  '200' ]= 0.019584 
    d['astral-true']['02'][  '400' ]= 0.013326 
    d['astral-true']['02'][ '1000' ]= 0.008744 
    d['astral-true']['05'][   '10' ]= 0.102914 
    d['astral-true']['05'][   '25' ]= 0.050004 
    d['astral-true']['05'][   '50' ]= 0.034578 
    d['astral-true']['05'][  '100' ]=  0.02292 
    d['astral-true']['05'][  '200' ]=  0.01416 
    d['astral-true']['05'][  '400' ]= 0.006658 
    d['astral-true']['05'][ '1000' ]= 0.003746 
    d['astral-true']['2'][   '10' ]= 0.099164 
    d['astral-true']['2'][   '25' ]= 0.057922 
    d['astral-true']['2'][   '50' ]= 0.038752 
    d['astral-true']['2'][  '100' ]= 0.025822 
    d['astral-true']['2'][  '200' ]= 0.019158 
    d['astral-true']['2'][  '400' ]= 0.011662 
    d['astral-true']['2'][ '1000' ]= 0.007492

    d['astral-true']['20'][   '10' ]= 0.120414
    d['astral-true']['20'][   '25' ]= 0.070836
    d['astral-true']['20'][   '50' ]= 0.038326
    d['astral-true']['20'][  '100' ]= 0.032906
    d['astral-true']['20'][  '200' ]= 0.017911
    d['astral-true']['20'][  '400' ]= 0.014992
    d['astral-true']['20'][ '1000' ]= 0.007494
    d['astral-true']['50'][   '10' ]= 0.235836
    d['astral-true']['50'][   '25' ]= 0.11542
    d['astral-true']['50'][   '50' ]= 0.072914
    d['astral-true']['50'][  '100' ]= 0.047916
    d['astral-true']['50'][  '200' ]= 0.035414
    d['astral-true']['50'][  '400' ]= 0.020412
    d['astral-true']['50'][ '1000' ]= 0.01208



    return d


def parse_filename(folder, fname):
    parts = fname.split('_')
    value = open(folder + "/" + fname).read()    
    try:
        value = float(value)

    except ValueError:
        print "couldn't read", fname
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
    d = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:-1)))
    d_raw = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
    for i in runs:
        d_raw[i.experiment][i.hgt][str(i.ngenes)].append(i.value)
    for experiment in sorted(d_raw.keys()):
        print experiment
        d1 = d[experiment]
        d1_raw = d_raw[experiment]
        for hgt in sorted(d1_raw.keys()):
            d2 = d1[hgt]
            d2_raw = d1_raw[hgt]
            for ngenes in sorted([int(i) for i in d2_raw.keys()]):
                print experiment, hgt, ngenes
                runs_averaged.append( Run(experiment, hgt, 0, ngenes, '', np.average(d2_raw[str(ngenes)])) )
                d2[str(ngenes)] = np.average(d2_raw[str(ngenes)])
#    d = add_astral(d)
    

    return d, runs_averaged, d_raw, runs

def printtable(methods_subset, hgtrates, genecounts, runs, d, separator=' ', endcharacter = ''):
    print hgtrates
    print genecounts
    for m in methods_subset:
        print d[m]

    def valordash(d, m, hgtrate, ngenes):
        try: 
            return '{:.1f}'.format(d[m][str(hgtrate)][str(ngenes)] * 100)
        except :
            return "-"
    
    for hgtrate in hgtrates:
        print "HGT"+separator+ "ngenes", separator.join([i.replace("missingbranchrate","") for i in methods_subset]) + endcharacter
        for ngenes in sorted([int(i) for i in genecounts]):            
            try:
                print str(hgtrate) + separator + str(ngenes) + separator + separator.join([valordash(d, m, str(hgtrate), str(ngenes)) for m in methods_subset]) + endcharacter
            except ValueError:
                print "ERROR", [(100 * d[m][str(hgtrate)][str(ngenes)]) for m in methods_subset]
                pass
                
        print
        print


def printdifferences(methods_subset, hgtrates, genecounts, runs, d, separator=' ', opener='', closer=''):
    print hgtrates
    print genecounts

    wqmc_better = []
    astral_better = []
    both_equal = []

    def valordash(d, m, hgtrate, ngenes):
        try: 
            return '{:.1f}'.format(d[m][str(hgtrate)][str(ngenes)] * 100)
        except :
            return "-"

    def differences(l0, l1):
        gt = []
        lt = []
        eq = []
        for i in l0:
            j = [k for k in l1 if k.index == i.index]
            if not len(j):
                continue
            j = j[0]
            if float(i.value) > float(j.value):
                gt.append((i, j))
                astral_better.append((i, j))
            elif float(i.value) < float(j.value):
                lt.append((i, j))
                wqmc_better.append((i, j))
            elif float(i.value) == float(j.value):
                eq.append((i, j))
                both_equal.append((i, j))
        return gt, eq, lt
            
    
    for hgtrate in hgtrates:
        print opener+"HGT"+separator+ "ngenes"+separator+methods_subset[0].replace("quartetscore","") +separator+ '=' + separator + methods_subset[1].replace("quartetscore","")+closer
        for ngenes in sorted([int(i) for i in genecounts]):            
            try:
                l0 = [i for i in runs if 
                      str(i.hgt) == str(hgtrate) and
                      str(i.experiment) == methods_subset[0] and
                      str(i.ngenes) == str(ngenes)]
                l1 = [i for i in runs if 
                      str(i.hgt) == str(hgtrate) and
                      str(i.experiment) == methods_subset[1] and
                      str(i.ngenes) == str(ngenes)]
                
                diff = differences(l0, l1)                                      
                print opener+str(hgtrate) + separator + str(ngenes) +separator + separator.join([str(len(i)) for i in diff])+closer
            except ValueError:
                print "ERROR", [(100 * d[m][str(hgtrate)][str(ngenes)]) for m in methods_subset]
                pass                
        print
        print
        
    gtc = 0
    ltc = 0
    eqc = 0
        
    print "|scenario| # times astral is better on MBR | # times astral = wQMC on MBR | # times wQMC is better on MBR "
    for i,j in wqmc_better:
        astralmbr = float(open('/home/vachasp2/scratch/branchrates/' + i.filename.replace('quartetscore', 'missingbranchrate')).read().strip())
        wqmcmbr = float(open('/home/vachasp2/scratch/branchrates/' + j.filename.replace('quartetscore', 'missingbranchrate')).read().strip())
#        print i.filename, j.filename
        if astralmbr < wqmcmbr:
            gtc+=1
        elif wqmcmbr < astralmbr :
            ltc += 1
        else:
            eqc += 1        

    # lower mbr is better

    print '|wQMC better quartet score |', gtc, " | ", eqc, " | ", ltc

    gtc = 0
    ltc = 0
    eqc = 0
        
    for i,j in astral_better:
        astralmbr = float(open('/home/vachasp2/scratch/branchrates/' + i.filename.replace('quartetscore', 'missingbranchrate')).read().strip())
        wqmcmbr = float(open('/home/vachasp2/scratch/branchrates/' + j.filename.replace('quartetscore', 'missingbranchrate')).read().strip())
#        print i.filename, j.filename
        if astralmbr < wqmcmbr:
            gtc+=1
        elif wqmcmbr < astralmbr :
            ltc += 1
        else:
            eqc += 1        
    print ""
    print "| ASTRAL better quartet score |", gtc, " | ", eqc, " | ", ltc

    gtc = 0
    ltc = 0
    eqc = 0
        
    for i,j in both_equal:
        astralmbr = float(open('/home/vachasp2/scratch/branchrates/' + i.filename.replace('quartetscore', 'missingbranchrate')).read().strip())
        wqmcmbr = float(open('/home/vachasp2/scratch/branchrates/' + j.filename.replace('quartetscore', 'missingbranchrate')).read().strip())
#        print i.filename, j.filename
        if astralmbr < wqmcmbr:
            gtc+=1
        elif wqmcmbr < astralmbr :
            ltc += 1
        else:
            eqc += 1        
    print "| Equal quartet scores |", gtc, " | ", eqc, " | ", ltc




def analyze(folder):
    d, runs, draw, runsraw = sort_files(folder)
    for experiment in sorted(d.keys()):
        print experiment
        d1 = d[experiment]
        for hgt in sorted(d1.keys()):
            d2 = d1[hgt]
            for ngenes in sorted([int(i) for i in d2.keys()]):
                print "|" + hgt + "\t|" + str(ngenes) + "\t|" + str(np.average(d2[str(ngenes)])) + "\t|" + str(np.std(d2[str(ngenes)])) + "\t|" + str(len(d2[str(ngenes)])) + "\t|"
                pass
        print 
    hgtrates = set([i.hgt for i in runs])
    methods = sorted(list(set([i.experiment for i in runs])))
    methods_true = sorted(list(set([i for i in methods if "true" in i])))
    methods_estimated = sorted(list(set([i for i in methods if "true" not in i])))
    methods_inv_true = sorted(list(set([i for i in methods_true if "invariants" in i])))
#    methods_inv_true = [methods_inv_true[-1]] + methods_inv_true[:-1]
#    methods_inv_estimated = sorted(list(set([i for i in methods_estimated if "invariants" in i])))

#    methods_sparsesample_estimated = sorted(list(set([i for i in methods_estimated if "sparsesample" in i and  i[-1] =='p']))) + [i for i in methods if 'wqmc-estimated' in i]
#    methods_sparsesample_true = sorted(list(set([i for i in methods_true if "sparsesample" in i and i[-1] =='p' ]))) + [i for i in methods if 'wqmc-true' in i]

    methods_important = ['astral', 'missingbranchratewqmc-estimated', 'missingbranchratewqmc-invariants',  'missingbranchratenjst-estimated', 'missingbranchrateastral-with-st-estimated',  'missingbranchrateastral-with-wqmc-estimated', 'astral-true', 'missingbranchratewqmc-true', 'missingbranchratewqmc-invariants-true', 'missingbranchratenjst-true','missingbranchrateastral-with-st-true',  'missingbranchrateastral-with-wqmc-true',]


    methods_quartetscores = ['quartetscoreastral-estimated', 'quartetscorewqmc-estimated', 'quartetscoreastral-true', 'quartetscorewqmc-true',]
    methods_quartetscores_estimated = ['quartetscoreastral-estimated', 'quartetscorewqmc-estimated']
    methods_quartetscores_true = ['quartetscoreastral-true', 'quartetscorewqmc-true']


    genecounts = set([i.ngenes for i in runs])


#    printtable(methods_estimated, hgtrates, genecounts, runs)
#    printtable(methods_true, hgtrates, genecounts, runs)

#    printtable(methods_inv_estimated, hgtrates, genecounts, runs)
#    printtable(methods_inv_true, hgtrates, genecounts, runs)

#    printtable(methods_sparsesample_estimated, hgtrates, genecounts, runs)
#    printtable(methods_sparsesample_true, hgtrates, genecounts, runs)

#    printtable(methods_quartetscores, hgtrates, genecounts, runs, d)

#    printdifferences(methods_quartetscores_true, hgtrates, genecounts, runsraw, draw, '|', '|', '|')
#    print 
#    print
#    printdifferences(methods_quartetscores_estimated, hgtrates, genecounts, runsraw, draw, '|', '|', '|')
#    print 
#    print    

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "python collate_branchrates.py ~/path/to/branchrates"
    folder = sys.argv[1]
    analyze(folder)
