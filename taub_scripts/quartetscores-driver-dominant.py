import sys
import copy
import quartetfile_dominantscores

infile = str(sys.argv[1])

outfile = str(sys.argv[2])

A = quartetfile_dominantscores.QuartetsInfo(infile,outfile)

l = len(A.quartet_dict())

print "read " + str(l) + " quartets from file " + infile

A.write_file()

print "new quartet invariant scores written to file " + outfile
