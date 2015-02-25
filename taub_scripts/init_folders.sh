# you should have in your basedir:
#  - whatever methods you're using; e.g. the ASTRAL folder, the wQMC folder, etc.
#  - the hgt-data folder (and whatever other data sets you're using)
#
# in principle you can put these wherever, but then you'll have to change some
# things in generate_qsub_script.py

mkdir ~/scratch/trees
mkdir ~/scratch/branchrates
mkdir ~/scratch/quartetscores
mkdir ~/scratch/timings

#make quartets script
cd quartets
cd thesis
make QUARTETS=1 CONTRACT_NUM=20000
cd ..
cp thesis/quart_bin triplets.soda2103
