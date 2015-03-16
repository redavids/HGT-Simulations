#import dendropy
#import matrixmaker
#import InvariantScores
import time
import sys

bigfile = str(sys.argv[1])

with open(bigfile) as f:
    content = f.readlines()
    f.close()

for j in range(len(content)):
    smallfilename = bigfile+str(j+1)+'.trees'
    g = open(smallfilename,'w')
    g.write(content[j])
    g.close()

         
