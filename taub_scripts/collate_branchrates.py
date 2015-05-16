from collections import defaultdict, namedtuple
import numpy as np
import os, sys

Run = namedtuple('Run', ['experiment', 'hgt', 'index', 'ngenes', 'filename', 'branchrate', 'quartetscore'])

def parse_filename(branchratefolder, quartetscorefolder, fname):
    parts = fname.split('_')
    branchrate = open(branchratefolder + "/" + fname).read()
    quartetscore = open(quartetscorefolder + "/" + fname.replace('missingbranchrate', 'quartetscore')).read()
    try:
        branchrate = float(branchrate)

    except ValueError:
        print "couldn't read", fname
        branchrate = -1
    try:
        quartetscore = float(quartetscore)
    except ValueError:
        print "couldn't read", fname.replace('missingbranchrate', 'quartetscore')
        branchrate = -1

#    parts[3] = int(parts[3])
#    if parts[3] == '.':
#        parts[3] = 1000
#    else:
#        parts[3] = int(parts[3].split('.')[1])
    
    return Run(parts[0], parts[1], parts[2], parts[3], fname, branchrate, quartetscore)
    
def sort_files(folders):
    runs = []
    for folder in folders:
        runs.extend([parse_filename(folder + '/branchrates/', folder + '/quartetscores', i) for i in os.listdir(folder + '/branchrates')])
    runs = [i for i in runs if i.branchrate >= 0]
    runs_averaged = []
    d = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:-1)))
    d_raw = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
    qd = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:-1)))
    qd_raw = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
    for i in runs:
        d_raw[i.experiment][i.hgt][str(i.ngenes)].append(i.branchrate)
        qd_raw[i.experiment][i.hgt][str(i.ngenes)].append(i.quartetscore)
        
    for experiment in sorted(d_raw.keys()):
        print experiment
        d1 = d[experiment]
        d1_raw = d_raw[experiment]
        qd1 = qd[experiment]
        qd1_raw = qd_raw[experiment]
        for hgt in sorted(d1_raw.keys()):
            d2 = d1[hgt]
            d2_raw = d1_raw[hgt]
            qd2 = qd1[hgt]
            qd2_raw = qd1_raw[hgt]
            for ngenes in sorted([int(i) for i in d2_raw.keys()]):
#                print experiment, hgt, ngenes
                runs_averaged.append( Run(experiment, hgt, 0, ngenes, '', np.average(d2_raw[str(ngenes)]), np.average(qd2_raw[str(ngenes)]) ) )
                d2[str(ngenes)] = (np.average(d2_raw[str(ngenes)]), np.std(d2_raw[str(ngenes)]))
                qd2[str(ngenes)] = np.average(qd2_raw[str(ngenes)])
#    d = add_astral(d)
    

    return d, qd, runs_averaged, d_raw, qd_raw, runs

def printgnuplot(methods_subset, hgtrates, genecounts, runs, d, qd, draw, qdraw, separator=' ', endcharacter = '', header = '', footer='', max_title=50, print_quartet_scores=False): 
    print hgtrates
    print genecounts
    for m in methods_subset:
        print d[m]
    color_names = sorted(list(set([i.split('.')[0] for i in methods_subset])))
    print color_names
    colors = dict(zip(color_names, range(len(color_names))))
    def valordash(d, m, hgtrate, ngenes):
        try: 
            if d[m][str(hgtrate)][str(ngenes)][0] * 100 == -100:
                return '-','-'
            return '{:.1f}'.format(d[m][str(hgtrate)][str(ngenes)][0] * 100), '{:.1f}'.format(d[m][str(hgtrate)][str(ngenes)][1] * 100)
        except :
            return "-"
    print methods_subset
    for hgtrate in hgtrates:
        print header
        for ngenes in sorted([int(i) for i in genecounts]):            
            methods_sub_subset = [m for m in methods_subset if d[m][str(hgtrate)][str(ngenes)] >= 0]
            if len(methods_sub_subset) == 0:
                continue
            print str(ngenes) + " " + hgtrate + " " + methods_sub_subset[0].split('.')[1]
            for m in methods_subset:
                if valordash(d, m, str(hgtrate), str(ngenes))[0] != '-':
                    print m.split('.')[0].replace('missingbranchrate', '') + separator + separator.join(valordash(d, m, str(hgtrate), str(ngenes))) + separator + str(colors[m.split('.')[0]])
#                    print m.split('.')[0].replace('missingbranchrate', '') + separator + separator.join([str(i) for i in draw[m][str(hgtrate)][str(ngenes)]]) + separator + str(colors[m.split('.')[0]])
            print
            print
             
        print footer
        print
        print



def printcsv(methods_subset, hgtrates, genecounts, runs, runsraw, d, qd, draw, qdraw, ofile):

    ofile.write(','.join(["experiment", "condition", "replicate", "ngenes", "branchrate", "quartetscore"]))
    ofile.write('\n')

    for run in runsraw:
        ofile.write(','.join([str(i) for i in [run.experiment, run.hgt, run.index, run.ngenes, run.branchrate, run.quartetscore]])) 
        ofile.write('\n')
 

def printtable(methods_subset, hgtrates, genecounts, runs, d, qd, separator=' ', endcharacter = '', header = '', footer='', max_title=50, print_quartet_scores=False): 
    print hgtrates
    print genecounts
    for m in methods_subset:
        print d[m]

    def valordash(d, m, hgtrate, ngenes):
        try: 
            return '{:.1f}'.format(d[m][str(hgtrate)][str(ngenes)] * 100)
        except :
            return "-"
    print methods_subset
    for hgtrate in hgtrates:
        print header
        print 'condition' + separator + separator.join([m.replace("missingbranchrate","") for m in methods_subset])
#        print "HGT"+separator+ "ngenes" + separator+ separator.join([i.replace("missingbranchrate","")[:max_title] for i in methods_subset]) + endcharacter
        for ngenes in sorted([int(i) for i in genecounts]):            
            try:
                print str(hgtrate) + separator + str(ngenes) + separator + separator.join([valordash(d, m, str(hgtrate), str(ngenes)) for m in methods_subset]) + endcharacter
            except ValueError:
                print "ERROR", [(100 * d[m][str(hgtrate)][str(ngenes)]) for m in methods_subset]
                pass
             
        print footer
        print
        print

def printtable_T(methods_subset, hgtrates, genecounts, runs, d, qd, separator=' ', endcharacter = '', header = '', footer='', max_title=50,  print_quartet_scores=False): 
    print hgtrates
    print genecounts
    for m in methods_subset:
        print d[m]

    def valordash(d, m, hgtrate, ngenes, ndec=1):
        try: 
            return ('{:.'+str(ndec)+'f}').format(d[m][str(hgtrate)][str(ngenes)][0] * 100)
        except :
            return "-"
    

    for hgtrate in hgtrates:

        print "\subsubsection{Branch Rate w/model condition " + hgtrate + "}"
        print header
#        print "HGT" + separator + separator.join([hgtrate] *  len(genecounts) ) + endcharacter
        print "ngenes" + separator + separator.join( [str(i) for i in sorted([int(i) for i in genecounts])]) + endcharacter

        for method in methods_subset:
            print method.replace("missingbranchrate","")[:max_title] + separator + separator.join([valordash(d, method, str(hgtrate), str(ngenes)) for ngenes in sorted([int(i) for i in genecounts])]) + endcharacter
             
        print footer
        print
        print

        if print_quartet_scores:

            print "\subsubsection{Quartet Score w/model condition " + hgtrate + "}"
            print header
        #        print "HGT" + separator + separator.join([hgtrate] *  len(genecounts) ) + endcharacter
            print "ngenes" + separator + separator.join( [str(i) for i in sorted([int(i) for i in genecounts])]) + endcharacter
            
            for method in methods_subset:
                print method.replace("missingbranchrate","")[:max_title] + separator + separator.join([valordash(qd, method, str(hgtrate), str(ngenes), 5) for ngenes in sorted([int(i) for i in genecounts])]) + endcharacter

            print footer





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
#        print opener+"HGT"+separator+ "ngenes"+separator+methods_subset[0].replace("quartetscore","") +separator+ '=' + separator + methods_subset[1].replace("quartetscore","")+closer
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



import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random

def plotscatter(methods, conditions, genecounts, runs, d, qd, draw, qdraw):
    pp = PdfPages('/home/vachasp2/plots.pdf')
    plotn = 0
#    plottypes = [321,322,323,324,325,326]
    plottypes = range(511, 516)
#    plottypes = [221, 224]
    fig = plt.figure(figsize=(8.5, 11), dpi=180)
    for condition in conditions:
        for ngenes in sorted([int(i) for i in genecounts]):
            xs = []
            ys = []
            sxs = []
            sys = []
            i = 0
            print condition, ngenes
            for m in methods:
                print m
                ys.append(draw[m][str(condition)][str(ngenes)])
                xs.append([i for k in ys[-1]])
                sys.extend(draw[m][str(condition)][str(ngenes)])
                sxs.extend([i + 0.1*(random.random() - 0.5) for k in ys[-1]])
                i += 1
            if len(ys) == 0:
                print "skipping"
                continue

            print sxs
            print sys
            xs = [i for i in xs if i]
            ys = [i for i in ys if i]
            print xs
            print ys
            thismethods=[j for j in methods if draw[j][str(condition)][str(ngenes)]]
#            plt.subplot2grid(gridsize, (plotn, 0), rowspan=1, colspan=1)
            plt.subplot(plottypes[plotn])
            plt.boxplot(ys, showmeans=True, labels=[i.replace("missingbranchrate", "").split('.')[0] for i in thismethods], positions=range(len(thismethods)), vert=False)
            plt.scatter(sys, sxs, marker='x')
            plt.title(thismethods[0].split('.')[1] + " " + condition + " " + str(ngenes))
            print "plotted", plottypes[plotn]
#            plt.yticks(rotation=45)
            plotn += 1
            if plotn==len(plottypes):
                plt.tight_layout()
                plotn = 0
                pp.savefig(fig)
                plt.close()
                fig = plt.figure(figsize=(8.5, 11), dpi=180)

    pp.close()

def analyze(folders, tpe, ofile):
    d, qd, runs, draw, qdraw, runsraw = sort_files(folders)
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

    if tpe == 'latex':
        printtable_T(methods, hgtrates, genecounts, runs, d, qd, separator='&', endcharacter='\\\\', header='\\begin{tabular}' + '{' + ' | '.join(['c'] * (len(genecounts) + 1)) + ' }', footer='\\end{tabular}\n\\\\\\hrule',max_title=50, print_quartet_scores=False)
    elif tpe == 'gnuplot':
        printgnuplot(methods, hgtrates, genecounts, runs, d, qd, draw, qdraw, separator=',', max_title=50, print_quartet_scores=False)
    elif tpe == 'scatter':
        plotscatter(methods, hgtrates, genecounts, runs, d, qd, draw, qdraw)
    elif tpe == 'csv':
        printcsv(methods, hgtrates, genecounts, runs, runsraw, d, qd, draw, qdraw, ofile)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "python collate_branchrates.py ~/path/to/branchrates"
    tpe = sys.argv[1]
    if tpe in ['csv']:
        folders=sys.argv[3:]
        ofile = sys.argv[2]
        analyze(folders, tpe, open(ofile, 'w'))
    else:
        folders = sys.argv[2:]
        analyze(folders, tpe, None)
