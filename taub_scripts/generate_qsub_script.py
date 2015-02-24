import sys
import itertools
import argparse
import os

basedir='/home/redavid2/phylogenetics'
scratchdir='/home/redavid2/scratch'

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
#PBS -M redavids@ncsu.edu
# By default, PBS scripts execute in your home directory, not the 
# directory from which they were submitted. The following line 
# places you in the directory from which the job was submitted.  
# run the program
#PBS -t 1-{njobs}

#env var for pbs_arrayid is 1 the first time 2 second, total number of tasks divided by tasks per ex reads first 30=tasks/job variable  inputs of param file

cd {basedir}

startindex=`expr "{tasksperjob}" "*" "${{PBS_ARRAYID}}"`
endindex=`expr "$startindex" "+" "{tasksperjob}" "-" "1"`

for ind in `seq $startindex $endindex`
do

#secs is variable in every bash script, double quotes {{ means is bash variable, won't string format it

starttime=$SECONDS

parameters=`cat {paramfile} | tail -n +${{ind}} | head -1`
parameterArray=($parameters)

identifier={jobname}_${{parameterArray[2]}}_${{parameterArray[0]}}_${{parameterArray[3]}}

genetreefilename=${{parameterArray[1]}}/${{parameterArray[0]}}/${{parameterArray[4]}}

quartetfilename={scratchdir}/quartets$identifier

fixedfilename={scratchdir}/fixed$identifier

treefilename={scratchdir}/trees/trees$identifier
timefilename={scratchdir}/timings/timing$identifier

branchratefilename={scratchdir}/branchrates/missingbranchrate$identifier
quartetscorefilename={scratchdir}/quartetscores/quartetscore$identifier


genetreesubsetfilename={scratchdir}/genetreesubset$identifier

speciestreefilename=${{parameterArray[1]}}/${{parameterArray[0]}}/s_tree.trees

head $genetreefilename -n${{parameterArray[3]}} > $genetreesubsetfilename

{methodcore}

compareTrees/compareTrees.missingBranchRate $speciestreefilename $treefilename > $branchratefilename

module load java

java -jar ASTRAL/astral.4.7.6.jar -q $treefilename -i $genetreesubsetfilename 2>&1 | tail -n1 | cut -f5 -d' ' > $quartetscorefilename

echo `expr $SECONDS - $starttime` > $timefilename

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
#PBS -M redavids@ncsu.edu
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

#tasks per job when > 1 when you have more than 1000 jobs only

tasks_per_job=1

def gen_param_file(paramfile, params):
    nreps = params['nreps']
    dirpath_labels = params['dirpaths']
    ngenes = params['ngenes']
    genetreetype = params['genetreetype']
    lines = [" ".join(i) for i in itertools.product(nreps, dirpath_labels, ngenes, [genetreetype])]
    s = "\n".join(lines)
    open(paramfile, "w").write(s)
    print len(nreps) * len(dirpath_labels) * len(ngenes) 
    return len(nreps) * len(dirpath_labels) * len(ngenes) 
    

def gen_main_qsub(jobname, params, nparams, paramfile):
    params['scratchdir'] = scratchdir
    params['basedir'] = basedir
    params['methodparams'].update(params)
    params['jobname'] = jobname
    params['njobs'] = nparams
    params['tasksperjob'] = tasks_per_job
    params['paramfile'] = paramfile
    params['methodcore'] = params['methodcore'].format(**(params['methodparams']))
    return mainformatstring.format(**params)
    
def gen_analyze_qsub(jobname):
    return collateformatstring.format(**{'jobname':jobname,
                                         'basedir':basedir,
                                         'scratchdir':scratchdir})

#dirpaths start from current directory to the directory containing data for experiment, followed by a short name
#Ruth was testing this with an older experiment so she commented out as was not necessary to re run all of them 
dirpaths = ['data/model.50.2000000.0.000001.0 0'] #,    
                    #'hgt-data/model.50.2000000.0.000001.0.000000002 02',
                    #'hgt-data/model.50.2000000.0.000001.0.000000005 05',
                    #'hgt-data/model.50.2000000.0.000001.0.00000002 2',
                    #'hgt-data/model.50.2000000.0.000001.0.0000002 20',
                    #'hgt-data/model.50.2000000.0.000001.0.0000005 50']


#names must not have underscores!

configs = {
    'wqmc-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'methodcore':wqmc_core,
        'methodparams':{
            'quartetscountfile':'quartetswqmc-estimated',
            'quartetgenerator':'quartets/quartet-controller.sh'
        }
    },
    'wqmc-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'methodcore':wqmc_core,
        'methodparams':{
            'quartetscountfile':'quartetswqmc-true',
            'quartetgenerator':'quartets/quartet-controller.sh'
        }
    },

    'astral-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'methodcore':astral_core,
        'mtehodparams':{}
    },
    'astral-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'methodcore':astral_core,
        'methodparams':{}
    },


    'astral-with-st-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'methodcore':astral_with_st_core,
        'methodparams':{}
    },
    'astral-with-st-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'methodcore':astral_with_st_core,
        'methodparams':{}
    },

    'astral-with-wqmc-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'methodcore':astral_with_wqmc_core,
        'methodparams':{
            'quartetgenerator':'bash ./astral_wqmc_qgen.sh',
            'quartetscountfile':'quartetswqmc-estimated'
        }
    },

    'astral-with-wqmc-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'methodcore':astral_with_wqmc_core,
        'methodparams': {
            'quartetscountfile':'quartetswqmc-true',
            'quartetgenerator':'bash ./astral_wqmc_qgen.sh'
        }
    },

    'njst-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'methodcore':njst_core,
        'methodparams':{}
    },
    'njst-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'methodcore':njst_core,
        'methodparams':{}
    },


    'wqmc-invariants': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'methodcore':wqmc_core,
        'methodparams':{
            'quartetscountfile':'quartetswqmc-estimated',
            'quartetgenerator':'bash ./invariant_qgen.sh 1',
        }
    },
    'wqmc-invariants-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'methodcore':wqmc_core,
        'methodparams':{
            'quartetscountfile':'quartetswqmc-true',
            'quartetgenerator':'bash ./invariant_qgen.sh 1'
        }
    },

    'wqmc-dominant': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'methodcore':wqmc_core,
        'methodparams':{
            'quartetscountfile':'quartetswqmc-estimated',
            'quartetgenerator':'bash ./dominant_qgen.sh',
        }
    },
    'wqmc-dominant-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'methodcore':wqmc_core,
        'methodparams': {
            'quartetscountfile':'quartetswqmc-true',
            'quartetgenerator':'bash ./dominant_qgen.sh',
        }
    },

}

for c in [0.1, 0.5, 2, 4, 10]:
    conf = configs['wqmc-invariants'].copy()
    conf['quartetgenerator'] = 'bash ./invariant_qgen.sh ' + str(c)
    configs['wqmc-invariants-c'+str(c)] = conf
    conf = configs['wqmc-invariants-true'].copy()
    conf['quartetgenerator'] = 'bash ./invariant_qgen.sh ' + str(c)
    configs['wqmc-invariants-c'+str(c)+'-true'] = conf


for p in [0.01, 0.1, 0.25, 0.5]:
    conf = configs['wqmc-estimated'].copy()
    conf['quartetgenerator'] = 'bash ./sparse_sample_qgen.sh ' + str(p)
    configs['wqmc-sparsesample-'+str(p)+'p'] = conf
    conf = configs['wqmc-true'].copy()
    conf['quartetgenerator'] = 'bash ./sparse_sample_qgen.sh ' + str(p)
    configs['wqmc-sparsesample-true-'+str(p)+'p'] = conf


if __name__ == "__main__":
    jobname = sys.argv[1]
    config = configs[jobname]
    paramfile = jobname + ".params"
    nparams = gen_param_file(paramfile, config)
    open(jobname + '.qsub', "w").write(gen_main_qsub(jobname, config, nparams, paramfile))
    open(jobname + '_analyze.qsub', "w").write(gen_analyze_qsub(jobname))
    if "--noqsub" not in sys.argv:
        os.system('qsub ' + jobname + '.qsub')
        os.system('qsub ' + jobname + '_analyze.qsub')
