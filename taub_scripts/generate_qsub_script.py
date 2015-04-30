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
#PBS -l cput=08:00:00  
# request 4 hours and 30 minutes of cpu time
#PBS -l walltime=08:00:00  
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m a
#PBS -q cse
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

starttime=$SECONDS

parameters=`cat {paramfile} | tail -n +${{ind}} | head -1`
parameterArray=($parameters)

replicate=${{parameterArray[0]}}
datafolder={datafolder}/${{parameterArray[1]}}
shortname=${{parameterArray[2]}}
ngenes=${{parameterArray[3]}}
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

done

exit 0
"""

use_truetree_core = """
cat $speciestreefilename > $treefilename

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



extend_bipartitions_core = """

module load R/3.1.2

mkdir $outputfolder/fulltrees


fulltreesfilename=$outputfolder/fulltrees/fulltrees$identifier
partialtreesfilename=$outputfolder/fulltrees/partialtrees$identifier
fastafilename=$outputfolder/fulltrees/fasta$identifier.fasta

sequencefilename={seqdatafolder}/${{parameterArray[1]}}/${{replicate}}/{seqfile}

rm $fulltreesfilename
touch $fulltreesfilename

rm $partialtreesfilename

echo "{siblingpairingexe} $genetreesubsetfilename >> $fulltreesfilename"

{siblingpairingexe} $genetreesubsetfilename >> $fulltreesfilename

{quartetgenerator} $genetreesubsetfilename $quartetfilename {scratchdir}/{quartetscountfile}-{dataset}_${{shortname}}_${{replicate}}_${{ngenes}}

cat $quartetfilename | sed s/"(("//g | sed s/"),("/"|"/g | sed s/")); "/":"/g | sed '/|/!d' > $fixedfilename

{wqmcexe} qrtt=$fixedfilename weights=on otre=$partialtreesfilename

cat $partialtreesfilename >> $fulltreesfilename


python /home/vachasp2/HGT-Simulations/concatenate.py $sequencefilename <(seq 1 1000) $ngenes > $fastafilename

$HOME/bin/fasttree -nt -gtr -nopr $fastafilename >> $fulltreesfilename

#{njstexe} $genetreesubsetfilename $partialtreesfilename
#cat $partialtreesfilename >> $fulltreesfilename

module load java

{astralexe} -i $genetreesubsetfilename -o $treefilename -e $fulltreesfilename

"""


extend_bipartitions_truetree_core = """

module load R/3.1.2

mkdir $outputfolder/fulltrees


fulltreesfilename=$outputfolder/fulltrees/fulltrees$identifier
partialtreesfilename=$outputfolder/fulltrees/partialtrees$identifier


rm $fulltreesfilename
touch $fulltreesfilename

rm $partialtreesfilename


module load java

{astralexe} -i $genetreesubsetfilename -o $treefilename -e $speciestreefilename

"""

extend_bipartitions_truetree_p2_core = """

module load R/3.1.2

mkdir $outputfolder/fulltrees


fulltreesfilename=$outputfolder/fulltrees/fulltrees$identifier
partialtreesfilename=$outputfolder/fulltrees/partialtrees$identifier


rm $fulltreesfilename
touch $fulltreesfilename

rm $partialtreesfilename


module load java

{astralexe} -i $genetreesubsetfilename -p2 -o $treefilename -e $speciestreefilename

"""



extend_bipartitions_completion_core = """

module load R/3.1.2

mkdir $outputfolder/fulltrees


fulltreesfilename=$outputfolder/fulltrees/fulltrees$identifier
completedtreesfilename=$outputfolder/fulltrees/completedtrees$identifier
partialtreesfilename=$outputfolder/fulltrees/partialtrees$identifier
fastafilename=$outputfolder/fulltrees/fasta$identifier.fasta

sequencefilename={seqdatafolder}/${{parameterArray[1]}}/${{replicate}}/{seqfile}

rm $fulltreesfilename
touch $fulltreesfilename

rm $partialtreesfilename

cat $genetreesubsetfilename >> $fulltreesfilename

echo "{siblingpairingexe} $genetreesubsetfilename >> $fulltreesfilename"

{siblingpairingexe} $genetreesubsetfilename >> $fulltreesfilename

python /home/vachasp2/HGT-Simulations/concatenate.py $sequencefilename <(seq 1 1000) $ngenes > $fastafilename

$HOME/bin/fasttree -nt -gtr -nopr $fastafilename >> $fulltreesfilename

#{njstexe} $genetreesubsetfilename $partialtreesfilename
#cat $partialtreesfilename >> $fulltreesfilename

module load python/2.7.8

python tree_completer.py $fulltreesfilename $completedtreesfilename


module load java

{astralexe} -i $genetreesubsetfilename -o $treefilename -e $completedtreesfilename

"""



extend_bipartitions_bootstrapped_completion_core = """

module load R/3.1.2

mkdir $outputfolder/fulltrees
mkdir $outputfolder/completedtrees
mkdir $outputfolder/bootstraps

completedtreesfilename=$outputfolder/completedtrees/completedtrees$identifier
fulltreesfilename=$outputfolder/fulltrees/fulltrees$identifier
bootstrapfilename=$outputfolder/bootstraps/boostrap$identifier
partialtreesfilename=$outputfolder/fulltrees/partialtrees$identifier


rm $fulltreesfilename
touch $fulltreesfilename

cat $genetreesubsetfilename >> $fulltreesfilename

for i in `seq 0 {nbootstraps}`
do

rm $bootstrapfilename
touch $bootstrapfilename

cat $genetreesubsetfilename | shuf | head -n {bootstrapsize} > $bootstrapfilename

{siblingpairingexe} $genetreesubsetfilename >> $fulltreesfilename

python /home/vachasp2/HGT-Simulations/concatenate.py $sequencefilename <(for i in `seq 1 1000`; do echo $((RANDOM % 1000 + 1)); done ) $ngenes > $fastafilename

$HOME/bin/fasttree -nt -gtr -nopr $fastafilename >> $fulltreesfilename


done

module load python/2.7.8

python tree_completer.py $fulltreesfilename $completedtreesfilename

module load java

{astralexe} -i $genetreesubsetfilename -o $treefilename -e $completedtreesfilename

"""



astral_core = """

module load java

{astralexe} -i $genetreesubsetfilename -o $treefilename

"""

astral_p2_core = """

module load java

{astralexe} -i $genetreesubsetfilename -p2 -o $treefilename

"""

astral_exact_core = """

module load java

{astralexe} -i $genetreesubsetfilename -x -o $treefilename

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

tasks_per_job=20

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

dirpaths = ['model.50.2000000.0.000001.0 0',    
                    'model.50.2000000.0.000001.0.000000002 02',
                    'model.50.2000000.0.000001.0.000000005 05',
                    'model.50.2000000.0.000001.0.00000002 2',
                    'model.50.2000000.0.000001.0.0000002 20',
                    'model.50.2000000.0.000001.0.0000005 50']

dirpaths_small = ['model.50.2000000.0.000001.0 0',    
                    'model.50.2000000.0.000001.0.0000002 2']


smidgen_dirpaths = ['100-taxa/100 100',
                    '500-taxa/100 500',
                    '1000-taxa/100 1000']

smidgen_20_dirpaths = ['100-taxa/20 100',
                       '500-taxa/20 500',
                       '1000-taxa/20 1000']

smidgen_50_dirpaths = ['100-taxa/50 100',
                       '500-taxa/50 500',
                       '1000-taxa/50 1000']


smidgen_75_dirpaths = ['100-taxa/75 100',
                       '500-taxa/75 500',
                       '1000-taxa/75 1000']

bansal_dirpaths = ['noise0 0000', 'noise500 0500', 'noise1000 1000', 'noise1500 1500', 'noise2000 2000',
                   'noise2500 2500', 'noise3000 3500', 'noise4000 4000', 'noise4500 4500', 
                   'noise5000 5000', 'noise5500 5500', 'noise6000 6000']
   
astral_2_dirpaths = ['model.10.2000000.0.000001 10-2M-1e6', 'model.100.2000000.0.000001 100-2M-1e6',# 'model.1000.2000000.0.000001 1000-2M-1e6',
                     'model.200.10000000.0.0000001 200-10M-1e7', 'model.200.10000000.0.000001 200-10M-1e6', 'model.200.2000000.0.0000001 200-2M-1e7', 
                     'model.200.2000000.0.000001 200-2M-1e6', 'model.200.500000.0.0000001 200-5M-1e7', 'model.200.500000.0.000001 200-5M-1e6',
                     'model.50.2000000.0.000001 50-2M-1e6', 'model.500.2000000.0.000001 500-2M-1e6']

astral_2_dirpaths_quick = ['model.10.2000000.0.000001 10-2M-1e6', 'model.50.2000000.0.000001 50-2M-1e6', 'model.100.2000000.0.000001 100-2M-1e6',
                           'model.200.2000000.0.000001 200-5M-1e7']

astral_2_dirpaths_10 = ['model.10.2000000.0.000001 10-2M-1e6']
astral_2_dirpaths_2002M1e6 = ['model.200.2000000.0.000001 200-2M-1e6']


astral_2_folder = '/home/vachasp2/scratch/astral-2-data/estimated/'

datasets = {
    'astral-2-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':astral_2_dirpaths,
        'datafolder':astral_2_folder,
        'ngenes':['10', '25', '50', '100'],
        'genetreetype':'estimatedgenetre',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylip'
    },

    'astral-2-quick-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':astral_2_dirpaths_quick,
        'datafolder':astral_2_folder,
        'ngenes':['10', '25', '50', '100'],
        'genetreetype':'estimatedgenetre',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylip'
    },

    'astral-2-10-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':astral_2_dirpaths_10,
        'datafolder':astral_2_folder,
        'ngenes':['10', '25', '50', '100'],
        'genetreetype':'estimatedgenetre',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylip'
    },

    'astral-2-2002M1e6-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':astral_2_dirpaths_2002M1e6,
        'datafolder':astral_2_folder,
        'ngenes':['10', '25', '50', '100'],
        'genetreetype':'estimatedgenetre',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylip'
    },

    'hgtdata-estimated': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetre',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylip'
    },
    'hgtdata-true': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetrees',
        'speciestreetype':'s_tree.trees',
        },


    'hgtdata-estimated-m01': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetrem01',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylipm01'
    },
    'hgtdata-true-m01': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetreesm01',
        'speciestreetype':'s_tree.trees',
        },


    'hgtdata-estimated-m05': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetrem05',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylipm05'
    },
    'hgtdata-true-m05': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetreesm05',
        'speciestreetype':'s_tree.trees',
        },


    'hgtdata-estimated-m10': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetrem10',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylipm10',
    },

    'hgtdata-estimated-m10-small': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths_small,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetrem10',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylipm10',
    },

    'hgtdata-true-m10': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetreesm10',
        'speciestreetype':'s_tree.trees',
        },


    'hgtdata-estimated-m15': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetrem15',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylipm15',
    },
    'hgtdata-true-m15': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetreesm15',
        'speciestreetype':'s_tree.trees',
        },


    'hgtdata-estimated-m20': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetrem20',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylipm20',
    },
    'hgtdata-true-m20': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetreesm20',
        'speciestreetype':'s_tree.trees',
        },


    'hgtdata-estimated-m25': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'datafolder':'hgt-data/',
        'seqdatafolder':'/home/vachasp2/scratch/alignments/',
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'estimatedgenetrem25',
        'speciestreetype':'s_tree.trees',
        'seqfile':'all-genes.phylipm25',
    },
    'hgtdata-true-m25': {
        'nreps':["%02d" % i for i in range(1,51)],
        'dirpaths':dirpaths,
        'ngenes':['10', '25', '50', '100', '200', '400', '1000'],
        'genetreetype':'truegenetreesm25',
        'speciestreetype':'s_tree.trees',
        },

    'smidgen-100': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_dirpaths,
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['1000'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },
    'smidgen-100-quick': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_dirpaths[:1],
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['10', '25', '50', '100', '200'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },

    'smidgen-20': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_20_dirpaths,
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['1000'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },
    'smidgen-20-quick': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_20_dirpaths[:1],
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['10', '25', '50', '100', '200'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },
    'smidgen-50': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_50_dirpaths,
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['1000'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },
    'smidgen-50-quick': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_50_dirpaths[:1],
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['10', '25', '50', '100', '200'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },
    'smidgen-75': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_75_dirpaths,
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['1000'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },
    'smidgen-75-quick': {
        'nreps':["%d" % i for i in range(0,31)],
        'dirpaths':smidgen_75_dirpaths[:1],
        'datafolder':'/home/vachasp2/smidgen/',
        'ngenes':['10', '25', '50', '100', '200'],
        'genetreetype':'source_trees',
        'speciestreetype':'model_tree',
        },
    'bansal-hw100' : {
        'nreps':["%d" % i for i in range(1,51)],
        'dirpaths':bansal_dirpaths,
        'datafolder':'/home/vachasp2/hw100',
        'ngenes':['10', '25', '50', '100', '1000'],
        'genetreetype':'shuftruegenetrees',
        'speciestreetype':'s_tree.trees'
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
    'astral-p2' : {
        'core':astral_p2_core,
        'params':{
        }
    },
    'astral-exact' : {
        'core':astral_exact_core,
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
    },

    'bipartition-extension-astral': {
        'core':extend_bipartitions_core,
        'params':{
            'quartetgenerator':'quartets/quartet-controller.sh',
            'quartetscountfile':'quartetswqmc-estimated-m10',
        }
    },

    'bipartition-extension-truetree-astral': {
        'core':extend_bipartitions_truetree_core,
        'params':{
        }
    },

    'bipartition-extension-truetree-astral-p2': {
        'core':extend_bipartitions_truetree_p2_core,
        'params':{
        }
    },    
    'bipartition-extension-completion-astral': {
        'core':extend_bipartitions_completion_core,
        'params':{
            'quartetgenerator':'quartets/quartet-controller.sh',
        }
    },    
    'bipartition-extension-bootstrapped-completion-astral': {
        'core':extend_bipartitions_bootstrapped_completion_core,
        'params':{
            'nbootstraps':'10',
            'bootstrapsize':'25'
            }
        },
    'truetree' : {
        'core':use_truetree_core,
        'params':{}
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

for m in [1, 2, 3, 4, 5, 10, 20, 50, 100, 200, 500]:
    conf = datasets['astral-2-estimated'].copy()
    conf['genetreetype'] = 'estimatedgenetrem%d' % m
    datasets['astral-2-estimated-m%03d' % m] = conf

    conf = datasets['astral-2-quick-estimated'].copy()
    conf['genetreetype'] = 'estimatedgenetrem%d' % m
    datasets['astral-2-quick-estimated-m%03d' % m] = conf

    conf = datasets['astral-2-10-estimated'].copy()
    conf['genetreetype'] = 'estimatedgenetrem%d' % m
    datasets['astral-2-10-estimated-m%03d' % m] = conf

    conf = datasets['astral-2-2002M1e6-estimated'].copy()
    conf['genetreetype'] = 'estimatedgenetrem%d' % m
    datasets['astral-2-2002M1e6-estimated-m%03d' % m] = conf


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "python generate_qsub_script.py methodname dataset"
        print "Methods:"
        for i in methods:
            print i
        print 
        print "Datasets:"
        for i in datasets:
            print i
        exit()
    methodname = sys.argv[1]
    datasetname = sys.argv[2]
    jobname = methodname + '.' + datasetname
    if methodname not in methods or datasetname not in datasets:
        print "Methods:"
        for i in methods:
            print i
        print 
        print "Datasets:"
        for i in datasets:
            print i
        exit()
    method = methods[methodname]
    dataset = datasets[datasetname]
    paramfile = jobname + ".params"
    nparams = gen_param_file(paramfile, dataset)
    open(jobname + '.qsub', "w").write(gen_main_qsub(jobname, methodname, datasetname, method, dataset, nparams, paramfile))
    open(jobname + '_analyze.qsub', "w").write(gen_analyze_qsub(jobname))
    if nparams/tasks_per_job > 1000:
        print "You're asking to create more than 1000 jobs! Try reducing the number of tasks or increasing tasks_per_job."
        exit() 
    if "--noqsub" not in sys.argv:
        os.system('qsub ' + jobname + '.qsub')
        os.system('qsub ' + jobname + '_analyze.qsub')
