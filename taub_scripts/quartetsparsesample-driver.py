import sys
import copy
import quartetfile_sparsesample

c = float(sys.argv[1])

infile = str(sys.argv[2])

outfile = str(sys.argv[3])

A = quartetfile_sparsesample.QuartetsInfo(c,infile,outfile)

l = len(A.quartet_dict())

print "read " + str(l) + " quartets from file " + infile

A.write_file()

print "new quartet invariant scores written to file " + outfile
