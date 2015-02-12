import sys
import itertools
import argparse
import os

basedir='/home/vachasp2/phylogenetics'
scratchdir='/home/vachasp2/scratch'

mainformatstring="""
# declare a name for this job to be sample_job
#PBS -N {jobname}
# request 1 node
#PBS -l nodes=1:ppn=1:taub
# request 4 hours and 30 minutes of cpu time
#PBS -l cput=04:00:00  
# request 4 hours and 30 minutes of cpu time
#PBS -l walltime=04:00:00  
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m a
#PBS -q cse
#PBS -l naccesspolicy=singleuser
# specify your email address
#PBS -M me@pranj.al
# By default, PBS scripts execute in your home directory, not the 
# directory from which they were submitted. The following line 
# places you in the directory from which the job was submitted.  
# run the program
#PBS -t 1-{njobs}

cd {basedir}

startindex=`expr "{tasksperjob}" "*" "${{PBS_ARRAYID}}"`
endindex=`expr "$startindex" "+" "{tasksperjob}" "-" "1"`

for ind in `seq $startindex $endindex`
do

parameters=`cat {paramfile} | tail -n +${{ind}} | head -1`
parameterArray=($parameters)

identifier={jobname}_${{parameterArray[2]}}_${{parameterArray[0]}}_${{parameterArray[3]}}

genetreefilename=${{parameterArray[1]}}/${{parameterArray[0]}}/${{parameterArray[4]}}

quartetfilename={scratchdir}/quartets$identifier

fixedfilename={scratchdir}/fixed$identifier

treefilename={scratchdir}/trees/trees$identifier

branchratefilename={scratchdir}/branchrates/missingbranchrate$identifier


genetreesubsetfilename={scratchdir}/genetreesubset$identifier

speciestreefilename=${{parameterArray[1]}}/${{parameterArray[0]}}/s_tree.trees

head $genetreefilename -n${{parameterArray[3]}} > $genetreesubsetfilename

{quartetgenerator} $genetreesubsetfilename $quartetfilename {scratchdir}/{quartetscountfile}_${{parameterArray[2]}}_${{parameterArray[0]}}_${{parameterArray[3]}}

cat $quartetfilename | sed s/"(("//g | sed s/"),("/"|"/g | sed s/")); "/":"/g | sed '/|/!d' > $fixedfilename

{method} qrtt=$fixedfilename {methodparams} otre=$treefilename

compareTrees/compareTrees.missingBranchRate $speciestreefilename $treefilename > $branchratefilename

done

exit 0
"""


collateformatstring="""
# declare a name for this job to be sample_job
#PBS -N {jobname}_analyze
# request 1 node
#PBS -l nodes=1:ppn=1:taub
# request 4 hours and 30 minutes of cpu time
#PBS -l cput=02:00:00  
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m bea
# specify your email address
#PBS -M me@pranj.al
#PBS -W {jobname}
# By default, PBS scripts execute in your home directory, not the 
# directory from which they were submitted. The following line 
# places you in the directory from which the job was submitted.  
# run the program

cd {basedir}

module load python/2.7.8

python collate_branchrates.py {scratchdir}/branchrates/ > {jobname}.out

exit 0
"""

tasks_per_job=50

def gen_param_file(paramfile, nreps, dirpath_labels, ngenes, genetreetype):
    lines = [" ".join(i) for i in itertools.product(nreps, dirpath_labels, ngenes, [genetreetype])]
    s = "\n".join(lines)
    open(paramfile, "w").write(s)
    print len(nreps) * len(dirpath_labels) * len(ngenes) 
    return len(nreps) * len(dirpath_labels) * len(ngenes) 
    

def gen_main_qsub(jobname, paramfile, nparams, quartetgenerator, method, methodparams, quartetsfile):
    basedir = '/home/vachasp2/phylogenetics/'
    params = {'jobname':jobname, 
              'njobs':str(nparams/tasks_per_job),
              'tasksperjob':str(tasks_per_job),
              'nparams':str(nparams),
              'basedir':basedir,
              'scratchdir':scratchdir,
              'paramfile':paramfile,
              'quartetgenerator':quartetgenerator,
              'method':method,
              'methodparams':methodparams,
              'quartetscountfile':quartetsfile}
    print params
    return mainformatstring.format(**params)
    
def gen_analyze_qsub(jobname):
    return collateformatstring.format(**{'jobname':jobname,
                                         'basedir':basedir,
                                         'scratchdir':scratchdir})

dirpaths = ['hgt-data/model.50.2000000.0.000001.0 0',    
                    'hgt-data/model.50.2000000.0.000001.0.000000002 02',
                    'hgt-data/model.50.2000000.0.000001.0.000000005 05',
                    'hgt-data/model.50.2000000.0.000001.0.00000002 2',
                    'hgt-data/model.50.2000000.0.000001.0.00000002 20',
                    'hgt-data/model.50.2000000.0.000001.0.00000002 50']

#names must not have underscores!
configs = {
    'wqmc-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 1',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 1',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c01': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.1',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c01-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.1',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c05': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.5',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c05-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.5',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c2': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 2',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c2-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 2',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c4': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 4',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c4-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 4',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c10': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 10',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-invariants-c10-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 10',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-dominant': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./dominant_qgen.sh',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-dominant-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./dominant_qgen.sh',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},

    'wqmc-sparsesample-01p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.01',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-sparsesample-true-01p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.01',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},

    'wqmc-sparsesample-10p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.1',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-sparsesample-true-10p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.1',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},

    'wqmc-sparsesample-25p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.25',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-sparsesample-true-25p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.25',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},

    'wqmc-sparsesample-50p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.5',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
    'wqmc-sparsesample-true-50p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.5',
        'method':'./wQMC/max-cut-tree',
        'methodparams':'weights=on'},
}

if __name__ == "__main__":
    jobname = sys.argv[1]
    config = configs[jobname]
    paramfile = jobname + ".params"
    nparams = gen_param_file(paramfile, config['nreps'], config['dirpaths'], config['ngenes'], config['genetreetype'])
    open(jobname + '.qsub', "w").write(gen_main_qsub(jobname, paramfile, nparams, config['quartetgenerator'], config['method'], config['methodparams'], config['quartetsfile']))
    open(jobname + '_analyze.qsub', "w").write(gen_analyze_qsub(jobname))
    if "--noqsub" not in sys.argv:
        os.system('qsub ' + jobname + '.qsub')
        os.system('qsub ' + jobname + '_analyze.qsub')
