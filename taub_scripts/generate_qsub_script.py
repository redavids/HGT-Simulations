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
quartetscorefilename={scratchdir}/quartetscores/quartetscore$identifier


genetreesubsetfilename={scratchdir}/genetreesubset$identifier

speciestreefilename=${{parameterArray[1]}}/${{parameterArray[0]}}/s_tree.trees

head $genetreefilename -n${{parameterArray[3]}} > $genetreesubsetfilename

{methodcore}

compareTrees/compareTrees.missingBranchRate $speciestreefilename $treefilename > $branchratefilename

module load java

java -jar ASTRAL/astral.4.7.6.jar -q $treefilename -i $genetreesubsetfilename 2>&1 | tail -n1 | cut -f5 -d' ' > $quartetscorefilename

done

exit 0
"""

wqmc_core = """

{quartetgenerator} $genetreesubsetfilename $quartetfilename {scratchdir}/{quartetscountfile}_${{parameterArray[2]}}_${{parameterArray[0]}}_${{parameterArray[3]}}

cat $quartetfilename | sed s/"(("//g | sed s/"),("/"|"/g | sed s/")); "/":"/g | sed '/|/!d' > $fixedfilename

./wQMC/max-cut-tree qrtt=$fixedfilename weights=on otre=$treefilename

"""

njst_core = """

module load R/3.1.2

./njst-package/njst $genetreesubsetfilename $treefilename

"""

astral_core = """

module load java

java -jar ./ASTRAL/astral.4.7.6.jar -i $genetreesubsetfilename -o $treefilename

"""

astral_with_st_core = """

module load java

java -jar ./ASTRAL/astral.4.7.6.jar -i $genetreesubsetfilename -e $speciestreefilename -o $treefilename

"""

astral_with_wqmc_core = """

{quartetgenerator} $genetreesubsetfilename $quartetfilename {scratchdir}/{quartetscountfile}_${{parameterArray[2]}}_${{parameterArray[0]}}_${{parameterArray[3]}}

cat $quartetfilename | sed s/"(("//g | sed s/"),("/"|"/g | sed s/")); "/":"/g | sed '/|/!d' > $fixedfilename

./wQMC/max-cut-tree qrtt=$fixedfilename weights=on otre=$treefilename-initial

module load java

java -jar ./ASTRAL/astral.4.7.6.jar -i $genetreesubsetfilename -e $treefilename-initial -o $treefilename

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

tasks_per_job=30

def gen_param_file(paramfile, nreps, dirpath_labels, ngenes, genetreetype):
    lines = [" ".join(i) for i in itertools.product(nreps, dirpath_labels, ngenes, [genetreetype])]
    s = "\n".join(lines)
    open(paramfile, "w").write(s)
    print len(nreps) * len(dirpath_labels) * len(ngenes) 
    return len(nreps) * len(dirpath_labels) * len(ngenes) 
    

def gen_main_qsub(jobname, paramfile, nparams, quartetgenerator, quartetsfile, methodcore):
    params = {'jobname':jobname, 
              'njobs':str(nparams/tasks_per_job),
              'tasksperjob':str(tasks_per_job),
              'nparams':str(nparams),
              'basedir':basedir,
              'scratchdir':scratchdir,
              'paramfile':paramfile,
              'quartetgenerator':quartetgenerator,
              'quartetscountfile':quartetsfile}
    params['methodcore'] = methodcore.format(**params)
    return mainformatstring.format(**params)
    
def gen_analyze_qsub(jobname):
    return collateformatstring.format(**{'jobname':jobname,
                                         'basedir':basedir,
                                         'scratchdir':scratchdir})

dirpaths = ['hgt-data/model.50.2000000.0.000001.0 0',    
                    'hgt-data/model.50.2000000.0.000001.0.000000002 02',
                    'hgt-data/model.50.2000000.0.000001.0.000000005 05',
                    'hgt-data/model.50.2000000.0.000001.0.00000002 2',
                    'hgt-data/model.50.2000000.0.000001.0.0000002 20',
                    'hgt-data/model.50.2000000.0.000001.0.0000005 50']

#names must not have underscores!
configs = {
    'wqmc-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'methodcore':wqmc_core},
    'wqmc-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'methodcore':wqmc_core
    },

    'astral-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'methodcore':astral_core},
    'astral-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'methodcore':astral_core
    },


    'astral-with-st-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'methodcore':astral_with_st_core},
    'astral-with-st-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'quartets/quartet-controller.sh',
        'methodcore':astral_with_st_core
    },

    'astral-with-wqmc-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./astral_wqmc_qgen.sh',
        'methodcore':astral_with_wqmc_core},
    'astral-with-wqmc-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./astral_wqmc_qgen.sh',
        'methodcore':astral_with_wqmc_core
    },

    'njst-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':None,
        'quartetgenerator':None,
        'methodcore':njst_core},
    'njst-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':None,
        'quartetgenerator':None,
        'methodcore':njst_core},


    'wqmc-invariants': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 1',
        'methodcore':wqmc_core},
    'wqmc-invariants-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 1',
        'methodcore':wqmc_core},
    'wqmc-invariants-c01': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.1',
        'methodcore':wqmc_core},
    'wqmc-invariants-c01-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.1',
        'methodcore':wqmc_core},
    'wqmc-invariants-c05': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.5',
        'methodcore':wqmc_core},
    'wqmc-invariants-c05-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 0.5',
        'methodcore':wqmc_core},
    'wqmc-invariants-c2': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 2',
        'methodcore':wqmc_core},
    'wqmc-invariants-c2-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 2',
        'methodcore':wqmc_core},
    'wqmc-invariants-c4': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 4',
        'methodcore':wqmc_core},
    'wqmc-invariants-c4-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 4',
        'methodcore':wqmc_core},
    'wqmc-invariants-c10': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./invariant_qgen.sh 10',
        'methodcore':wqmc_core},
    'wqmc-invariants-c10-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./invariant_qgen.sh 10',
        'methodcore':wqmc_core},
    'wqmc-dominant': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./dominant_qgen.sh',
        'methodcore':wqmc_core},
    'wqmc-dominant-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./dominant_qgen.sh',
        'methodcore':wqmc_core},

    'wqmc-sparsesample-01p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.01',
        'methodcore':wqmc_core},
    'wqmc-sparsesample-true-01p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.01',
        'methodcore':wqmc_core},

    'wqmc-sparsesample-10p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.1',
        'methodcore':wqmc_core},
    'wqmc-sparsesample-true-10p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.1',
        'methodcore':wqmc_core},

    'wqmc-sparsesample-25p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.25',
        'methodcore':wqmc_core},
    'wqmc-sparsesample-true-25p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.25',
        'methodcore':wqmc_core},

    'wqmc-sparsesample-50p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'quartetsfile':'quartetswqmc-estimated',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.5',
        'methodcore':wqmc_core},
    'wqmc-sparsesample-true-50p': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'quartetsfile':'quartetswqmc-true',
        'quartetgenerator':'bash ./sparse_sample_qgen.sh 0.5',
        'methodcore':wqmc_core},
}

if __name__ == "__main__":
    jobname = sys.argv[1]
    config = configs[jobname]
    paramfile = jobname + ".params"
    nparams = gen_param_file(paramfile, config['nreps'], config['dirpaths'], config['ngenes'], config['genetreetype'])
    open(jobname + '.qsub', "w").write(gen_main_qsub(jobname, paramfile, nparams, config['quartetgenerator'], config['quartetsfile'], config['methodcore']))
    open(jobname + '_analyze.qsub', "w").write(gen_analyze_qsub(jobname))
    if "--noqsub" not in sys.argv:
        os.system('qsub ' + jobname + '.qsub')
        os.system('qsub ' + jobname + '_analyze.qsub')
