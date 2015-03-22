from collections import defaultdict, namedtuple
import numpy as np
import os, sys

Run = namedtuple('Run', ['experiment', 'hgt', 'index', 'ngenes', 'filename', 'value'])

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
    
def sort_files(folders):
    runs = []
    for folder in folders:
        runs.extend([parse_filename(folder, i) for i in os.listdir(folder)])
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
#                print experiment, hgt, ngenes
                runs_averaged.append( Run(experiment, hgt, 0, ngenes, '', np.average(d2_raw[str(ngenes)])) )
                d2[str(ngenes)] = np.average(d2_raw[str(ngenes)])
#    d = add_astral(d)
    

    return d, runs_averaged, d_raw, runs

def printtable(methods_subset, hgtrates, genecounts, runs, d, separator=' ', endcharacter = '', header = '', footer='', max_title=50): 
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
        print header
        print "HGT"+separator+ "ngenes" + separator+ separator.join([i.replace("missingbranchrate","")[:max_title] for i in methods_subset]) + endcharacter
        for ngenes in sorted([int(i) for i in genecounts]):            
            try:
                print str(hgtrate) + separator + str(ngenes) + separator + separator.join([valordash(d, m, str(hgtrate), str(ngenes)) for m in methods_subset]) + endcharacter
            except ValueError:
                print "ERROR", [(100 * d[m][str(hgtrate)][str(ngenes)]) for m in methods_subset]
                pass
             
        print footer
        print
        print

def printtable_T(methods_subset, hgtrates, genecounts, runs, d, separator=' ', endcharacter = '', header = '', footer='', max_title=50): 
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
        print header
        print "HGT" + separator + separator.join([hgtrate] *  len(genecounts) ) + endcharacter
        print "ngenes" + separator + separator.join( [str(i) for i in sorted([int(i) for i in genecounts])]) + endcharacter

        for method in methods_subset:
            print method.replace("missingbranchrate","")[:max_title] + separator + separator.join([valordash(d, method, str(hgtrate), str(ngenes)) for ngenes in sorted([int(i) for i in genecounts])]) + endcharacter
             
        print footer
        print
        print



def printraw(methods_subset, hgtrates, genecounts, runs, d, separator=',', endcharacter = ''): 
    print 'experiment,hgtrate,ngenes,index,missingbranchrate'
    for r in runs:
        print separator.join([r.experiment, r.hgt, r.ngenes, r.index, r.value])

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




def analyze(folders):
    d, runs, draw, runsraw = sort_files(folders)
    for experiment in sorted(d.keys()):
        print experiment
        d1 = d[experiment]
        d1raw = draw[experiment]
        for hgt in sorted(d1.keys()):
            d2 = d1[hgt]
            d2raw = d1raw[hgt]
            for ngenes in sorted([int(i) for i in d2.keys()]):
                print "|" + hgt + "\t|" + str(ngenes) + "\t|" + str(np.average(d2[str(ngenes)])) + "\t|" + str(np.std(d2[str(ngenes)])) + "\t|" + str(len(d2raw[str(ngenes)])) + "\t|"
                pass
        print 

    hgtrates = set([i.hgt for i in runs])
    genecounts = set([i.ngenes for i in runs])
    methods = sorted(list(set([i.experiment for i in runs])))
    methods_true = sorted(list(set([i for i in methods if "true" in i])))
    methods_estimated = sorted(list(set([i for i in methods if "true" not in i])))
    methods_inv_true = sorted(list(set([i for i in methods_true if "invariants" in i])))

    printtable_T(methods, hgtrates, genecounts, runs, d, separator='&', endcharacter='\\\\', header='\\begin{tabular}' + '{' + ' | '.join(['c'] * (len(genecounts) + 1)) + ' }', footer='\\end{tabular}\n\\\\\\hrule',max_title=5)

#    printraw(methods, hgtrates, genecounts, runsraw, draw, separator=',', endcharacter = '')

#    methods_inv_true = [methods_inv_true[-1]] + methods_inv_true[:-1]
#    methods_inv_estimated = sorted(list(set([i for i in methods_estimated if "invariants" in i])))

#    methods_sparsesample_estimated = sorted(list(set([i for i in methods_estimated if "sparsesample" in i and  i[-1] =='p']))) + [i for i in methods if 'wqmc-estimated' in i]
#    methods_sparsesample_true = sorted(list(set([i for i in methods_true if "sparsesample" in i and i[-1] =='p' ]))) + [i for i in methods if 'wqmc-true' in i]

    methods_important = ['astral', 'missingbranchratewqmc-estimated', 'missingbranchratewqmc-invariants',  'missingbranchratenjst-estimated', 'missingbranchrateastral-with-st-estimated',  'missingbranchrateastral-with-wqmc-estimated', 'astral-true', 'missingbranchratewqmc-true', 'missingbranchratewqmc-invariants-true', 'missingbranchratenjst-true','missingbranchrateastral-with-st-true',  'missingbranchrateastral-with-wqmc-true',]


    methods_quartetscores = ['quartetscoreastral-estimated', 'quartetscorewqmc-estimated', 'quartetscoreastral-true', 'quartetscorewqmc-true',]
    methods_quartetscores_estimated = ['quartetscoreastral-estimated', 'quartetscorewqmc-estimated']
    methods_quartetscores_true = ['quartetscoreastral-true', 'quartetscorewqmc-true']



#    printtable(methods_estimated, hgtrates, genecounts, runs)
#    printtable(methods_true, hgtrates, genecounts, runs)

#    printtable(methods_inv_estimated, hgtrates, genecounts, runs)
#    printtable(methods_inv_true, hgtrates, genecounts, runs)

#    printtable(methods_sparsesample_estimated, hgtrates, genecounts, runs)
#    printtable(methods_sparsesample_true, hgtrates, genecounts, runs)


#    printdifferences(methods_quartetscores_true, hgtrates, genecounts, runsraw, draw, '|', '|', '|')
#    print 
#    print
#    printdifferences(methods_quartetscores_estimated, hgtrates, genecounts, runsraw, draw, '|', '|', '|')
#    print 
#    print    

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "python collate_branchrates.py ~/path/to/branchrates"
    folders = sys.argv[1:]
    analyze(folders)
