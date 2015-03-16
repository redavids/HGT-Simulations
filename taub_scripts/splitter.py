#import dendropy
#import matrixmaker
#import InvariantScores
import time
import sys

#bigfile = str(sys.argv[1])

directory = str(sys.argv[1])

chompthefile = str(sys.argv[2])

with open(directory+chompthefile) as f:
    content = f.readlines()
    f.close()

for j in range(len(content)):
    smallfilename = directory + 'g_trees'+str(j+1)+'.trees'
    g = open(smallfilename,'w')
    g.write(content[j])
    g.close()

         
