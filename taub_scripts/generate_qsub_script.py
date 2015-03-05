import sys
import itertools
import argparse
import os
from settings import *

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
#PBS -M {email}
# By default, PBS scripts execute in your home directory, not the 
# directory from which they were submitted. The following line 
# places you in the directory from which the job was submitted.  
# run the program
#PBS -t 0-{njobs}

cd {basedir}

startindex=`expr "{tasksperjob}" "*" "${{PBS_ARRAYID}}"`
endindex=`expr "$startindex" "+" "{tasksperjob}" "-" "1"`

for ind in `seq $startindex $endindex`
do

starttime=$SECONDS

parameters=`cat {paramfile} | tail -n +${{ind}} | head -1`
parameterArray=($parameters)

replicate=${{parameterArray[0]}}
datafolder=${{parameterArray[1]}}
shortname=${{parameterArray[2]}}
ngenes=${{parameterArray[3]}}
genetree=${{parameterArray[4]}}
speciestree=${{parameterArray[5]}}

genetreefilename=${{datafolder}}/${{replicate}}/${{genetree}}
speciestreefilename=${{datafolder}}/${{replicate}}/${{speciestree}}

identifier={jobname}_${{shortname}}_${{replicate}}_${{ngenes}}
outputfolder={scratchdir}/{datasetname}/{methodname}/

mkdir -p $outputfolder
mkdir -p $outputfolder/quartets
mkdir -p $outputfolder/trees
mkdir -p $outputfolder/genetreesubsets

quartetfilename={scratchdir}/quartets$identifier

fixedfilename={scratchdir}/fixed$identifier

treefilename={scratchdir}/trees/trees$identifier
timefilename={scratchdir}/timings/timing$identifier

branchratefilename={scratchdir}/branchrates/missingbranchrate$identifier
quartetscorefilename={scratchdir}/quartetscores/quartetscore$identifier


genetreesubsetfilename={scratchdir}/genetreesubset$identifier

head $genetreefilename -n${{ngenes}} > $genetreesubsetfilename

{methodcore}

{comparer} $speciestreefilename $treefilename > $branchratefilename

module load java

{astralexe} -q $treefilename -i $genetreesubsetfilename 2>&1 | tail -n1 | cut -f5 -d' ' > $quartetscorefilename

echo `expr $SECONDS - $starttime` > $timefilename

done

exit 0
"""

wqmc_core = """

{quartetgenerator} $genetreesubsetfilename $quartetfilename {scratchdir}/{quartetscountfile}-{dataset}_${{shortname}}_${{replicate}}_${{ngenes}}

cat $quartetfilename | sed s/"(("//g | sed s/"),("/"|"/g | sed s/")); "/":"/g | sed '/|/!d' > $fixedfilename

{wqmcexe} qrtt=$fixedfilename weights=on otre=$treefilename

"""

njst_core = """

module load R/3.1.2

{njstexe} $genetreesubsetfilename $treefilename

"""

sibling_pairing_core = """

module load python/2.7.8

{siblingpairingexe} $genetreesubsetfilename > $treefilename

"""


astral_core = """

module load java

{astralexe} -i $genetreesubsetfilename -o $treefilename

"""

astral_with_st_core = """

module load java

{astralexe} -i $genetreesubsetfilename -e $speciestreefilename -o $treefilename

"""

astral_with_wqmc_core = """

{quartetgenerator} $genetreesubsetfilename $quartetfilename {scratchdir}/{quartetscountfile}-{dataset}_${{shortname}}_${{replicate}}_${{ngenes}}

cat $quartetfilename | sed s/"(("//g | sed s/"),("/"|"/g | sed s/")); "/":"/g | sed '/|/!d' > $fixedfilename

{wqmcexe} qrtt=$fixedfilename weights=on otre=$treefilename-initial

module load java

{astralexe} -i $genetreesubsetfilename -e $treefilename-initial -o $treefilename

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
#PBS -M {email}
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

def gen_param_file(paramfile, params):
    nreps = params['nreps']
    dirpath_labels = params['dirpaths']
    ngenes = params['ngenes']
    genetreetype = params['genetreetype']
    speciestreetype = params['speciestreetype']
    lines = [" ".join(i) for i in itertools.product(nreps, dirpath_labels, ngenes, [genetreetype], [speciestreetype])]
    s = "\n".join(lines)
    open(paramfile, "w").write(s)
    print len(nreps) * len(dirpath_labels) * len(ngenes) 
    return len(nreps) * len(dirpath_labels) * len(ngenes) 
    

def gen_main_qsub(jobname, methodname, datasetname, method, dataset, nparams, paramfile):
    params = {}
    params['jobname'] = jobname
    params['method'] = methodname
    params['dataset'] = datasetname
    params['njobs'] = nparams/tasks_per_job
    params['tasksperjob'] = tasks_per_job
    params['paramfile'] = paramfile
    
    params.update(globalparams)
    params.update(method)
    params.update(method['params'])
    params.update(dataset)
    print params
    params['methodcore'] = method['core'].format(**params)    

    return mainformatstring.format(**params)
    
def gen_analyze_qsub(jobname):
    params = {'jobname':jobname}
    params.update(globalparams)
    return collateformatstring.format(**params)

dirpaths = ['hgt-data/model.50.2000000.0.000001.0 0',    
                    'hgt-data/model.50.2000000.0.000001.0.000000002 02',
                    'hgt-data/model.50.2000000.0.000001.0.000000005 05',
                    'hgt-data/model.50.2000000.0.000001.0.00000002 2',
                    'hgt-data/model.50.2000000.0.000001.0.0000002 20',
                    'hgt-data/model.50.2000000.0.000001.0.0000005 50']


datasets = {
    'hgtdata-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'speciestreetype':'s_tree.trees',
    },
    'hgtdata-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'speciestreetype':'s_tree.trees',
    }
}

methods = {
    'wqmc' : {
        'core':wqmc_core,
        'params':{
            'quartetscountfile':'quartetswqmc',
            'quartetgenerator':'quartets/quartet-controller.sh',
        }
    },
    'siblingpairing' : {
        'core':sibling_pairing_core,
        'params':{
        }
    },
    'astral' : {
        'core':astral_core,
        'params':{
        }
    },
    'astral-with-st' : {
        'core':astral_with_st_core,
        'params':{}
    },
    'astral-with-wqmc' : {
        'core':astral_with_wqmc_core,
        'params':{
            'quartetgenerator':'bash ./astral_wqmc_qgen.sh',
            'quartetscountfile':'quartetswqmc-estimated'
        }
    },
    'njst': {
        'core':njst_core,
        'params':{}
    },
    'wqmc-invariants' : {
        'core':wqmc_core,
        'params':{
            'quartetscountfile':'quartetswqmc-estimated',
            'quartetgenerator':'bash ./invariant_qgen.sh 1',
        }
    },
    'wqmc-dominant': {
        'core':wqmc_core,
        'params':{
            'quartetscountfile':'quartetswqmc-estimated',
            'quartetgenerator':'bash ./dominant_qgen.sh',
        }
    }
    
}

for c in [0.1, 0.5, 2, 4, 10]:
    conf = methods['wqmc-invariants'].copy()
    conf['quartetgenerator'] = 'bash ./invariant_qgen.sh ' + str(c)
    methods['wqmc-invariants-c'+str(c)] = conf


for p in [0.01, 0.1, 0.25, 0.5]:
    conf = methods['wqmc'].copy()
    conf['quartetgenerator'] = 'bash ./sparse_sample_qgen.sh ' + str(p)
    methods['wqmc-sparsesample-'+str(p)+'p'] = conf

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "python generate_qsub_script.py methodname dataset"
    methodname = sys.argv[1]
    datasetname = sys.argv[2]
    jobname = methodname + '.' + datasetname
    method = methods[methodname]
    dataset = datasets[datasetname]
    paramfile = jobname + ".params"
    nparams = gen_param_file(paramfile, dataset)
    open(jobname + '.qsub', "w").write(gen_main_qsub(jobname, methodname, datasetname, method, dataset, nparams, paramfile))
    open(jobname + '_analyze.qsub', "w").write(gen_analyze_qsub(jobname))
    if "--noqsub" not in sys.argv:
        os.system('qsub ' + jobname + '.qsub')
        os.system('qsub ' + jobname + '_analyze.qsub')
