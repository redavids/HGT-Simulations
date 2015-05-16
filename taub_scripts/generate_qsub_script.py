import sys
import itertools
import argparse
import os
from settings import *

parallelrunner="""
# declare a name for this job to be sample_job
#PBS -N {jobname}
# request 1 node
#PBS -l nodes=1:ppn={parallelism}
# request 4 hours and 30 minutes of cpu time
#PBS -l cput=720:00:00  
# request 4 hours and 30 minutes of cpu time
#PBS -l walltime=36:00:00  
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m a
#PBS -q tallis
#PBS -l naccesspolicy=singleuser
# merge output and error files into output file
#PBS -j oe
# put output and error files in log directory
#PBS -o {scratchdir}/logs/
# specify your email address
#PBS -M {email}
# By default, PBS scripts execute in your home directory, not the 
# directory from which they were submitted. The following line 
# places you in the directory from which the job was submitted.  
# run the program
#PBS -t 0-{njobs}

source ~/.bashrc

cd {basedir}

startindex=`expr "{tasksperjob}" "*" "${{PBS_ARRAYID}}"`
endindex=`expr "$startindex" "+" "{tasksperjob}" "-" "1"`

for ind in `seq $startindex $endindex`
do

bash {jobname} $ind &

done

wait
exit 0

"""

mainformatstring="""
ind=$1

starttime=$SECONDS

parameters=`cat {paramfile} | tail -n +${{ind}} | head -1`
parameterArray=($parameters)

ngenes=${{parameterArray[0]}}
replicate=${{parameterArray[1]}}
datafolder={datafolder}/${{parameterArray[2]}}
shortname=${{parameterArray[3]}}
genetree=${{parameterArray[4]}}
speciestree=${{parameterArray[5]}}

genetreefilename=${{datafolder}}/${{replicate}}/${{genetree}}
speciestreefilename=${{datafolder}}/${{replicate}}/${{speciestree}}

identifier={jobname}_${{shortname}}_${{replicate}}_${{ngenes}}
outputfolder={scratchdir}/{dataset}/{method}/
resultsfolder={resultsdir}/{dataset}/{method}/

mkdir -p $outputfolder
mkdir -p $outputfolder/quartets
mkdir -p $outputfolder/trees
mkdir -p $outputfolder/branchrates
mkdir -p $outputfolder/quartetscores
mkdir -p $outputfolder/timings
mkdir -p $outputfolder/genetreesubsets


mkdir -p $resultsfolder
mkdir -p $resultsfolder/trees
mkdir -p $resultsfolder/timings
mkdir -p $resultsfolder/branchrates
mkdir -p $resultsfolder/quartetscores

quartetfilename={scratchdir}/quartets$identifier

fixedfilename={scratchdir}/fixed$identifier

treefilename=$resultsfolder/trees/trees$identifier
timefilename=$resultsfolder/timings/timing$identifier

branchratefilename=$resultsfolder/branchrates/missingbranchrate$identifier
quartetscorefilename=$resultsfolder/quartetscores/quartetscore$identifier


genetreesubsetfilename=$outputfolder/genetreesubsets/genetreesubset$identifier

cat $genetreefilename  | head -n${{ngenes}} > $genetreesubsetfilename

{methodcore}

{comparer} <(cat $speciestreefilename <(echo ';')) $treefilename > $branchratefilename

module load java

{astralexe} -q $treefilename -i $genetreesubsetfilename 2>&1 | tail -n1 | cut -f5 -d' ' > $quartetscorefilename

echo `expr $SECONDS - $starttime` > $timefilename

"""

collateformatstring="""
# declare a name for this job to be sample_job
#PBS -N {jobname}_analyze
# request 1 node
#PBS -l nodes=1:ppn=1
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



def gen_param_file(paramfile, params):
    nreps = params['nreps']
    dirpath_labels = params['dirpaths']
    ngenes = params['ngenes']
    genetreetype = params['genetreetype']
    speciestreetype = params['speciestreetype']
    lines = [" ".join(i) for i in itertools.product(ngenes, nreps, dirpath_labels, [genetreetype], [speciestreetype])]
    s = "\n".join(lines)
    open(paramfile, "w").write(s)
    print len(nreps) * len(dirpath_labels) * len(ngenes) 
    return len(nreps) * len(dirpath_labels) * len(ngenes) 
    

def gen_main_qsub(jobname, methodname, datasetname, method, dataset, nparams, paramfile, tasks_per_job):
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

def gen_parallel_qsub(jobname, methodname, datasetname, method, dataset, nparams, paramfile, tasks_per_job):
    params = {}
    params['jobname'] = jobname
    params['method'] = methodname
    params['dataset'] = datasetname
    params['njobs'] = nparams/tasks_per_job
    params['tasksperjob'] = tasks_per_job
    params['paramfile'] = paramfile
    params['parallelism'] = min(tasks_per_job, 20)

    params.update(globalparams)
    params.update(method)
    params.update(method['params'])
    params.update(dataset)
    print params
    params['methodcore'] = method['core'].format(**params)    

    return parallelrunner.format(**params)


    
def gen_analyze_qsub(jobname):
    params = {'jobname':jobname}
    params.update(globalparams)
    return collateformatstring.format(**params)

from methods import methods

from datasets import datasets

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run phylogenetics pipeline")
    parser.add_argument('-m', '--methodname', dest = 'methodnames', choices=sorted(methods.keys()), action='append', required=True)
    parser.add_argument('-d', '--datasetname', dest = 'datasetnames', choices=sorted(datasets.keys()), action='append', required=True)
    parser.add_argument('-q', '--qsub', dest='qsub', action='store_true')
    parser.add_argument('-t', '--tasksperjob', dest='tasks_per_job', type=int, default=10)


    args = parser.parse_args(sys.argv[1:])

    print args

    for methodname in args.methodnames:
        for datasetname in args.datasetnames:
            jobname = methodname + '.' + datasetname
            method = methods[methodname]
            dataset = datasets[datasetname]
            paramfile = jobname + ".params"
            nparams = gen_param_file(paramfile, dataset)
            open(jobname, "w").write(gen_main_qsub(jobname, methodname, datasetname, method, dataset, nparams, paramfile, args.tasks_per_job))
            open(jobname + '.qsub', "w").write(gen_parallel_qsub(jobname, methodname, datasetname, method, dataset, nparams, paramfile, args.tasks_per_job))

            open(jobname + '_analyze.qsub', "w").write(gen_analyze_qsub(jobname))
            if nparams/args.tasks_per_job > 1000:
                print "You're asking to create more than 1000 jobs! Try reducing the number of tasks or increasing tasks_per_job."
                exit() 
            if args.qsub:
                os.system('qsub ' + jobname + '.qsub')
#        os.system('qsub ' + jobname + '_analyze.qsub')
    
