import sys
from collections import defaultdict
import numpy as np

genecount=None
if len(sys.argv) > 3:
	genecount=int(sys.argv[3])

subset=None
if len(sys.argv) > 2:
	subset = set()
	gc=0
	for line in open(sys.argv[2]):
		subset.add(int(line.strip()))
		gc += 1
		if genecount is not None and gc >= genecount:
			break


s=0
i=0
genes = []

for line in open(sys.argv[1]):
	s = line.split()
        if not s[1].isdigit():
		genes[-1][s[0]] = s[1]
	else:
		i += 1
		genes.append({})
	

sys.stderr.write("subsetlength" +  str(len(subset)) + '\n')
sys.stderr.write("genelength" +  str(len(genes)) + '\n')

i=0
alg=defaultdict(lambda:[])
for s in subset:
	i += 1
	for x, g in genes[s].items():
		alg[x].append(g)

sys.stderr.write("i" +  str(i) + '\n')

l = 0
for k,vl in alg.iteritems():
	v=''.join(vl)
	if l == 0:
		l = len(v)
	if l != len(v):
		raise Exception("Key %s has %d columns instead of %d" %(k,len(v),l))
	print ">"+k
	print v
