import dendropy
import numpy as np

def delete_taxa(t, n):
    t.prune_taxa([i.taxon for i in np.random.choice(t.leaf_nodes(), size=n)])
    

import sys
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "delete_taxa.py input ndelete output"
        exit()
    inp = sys.argv[1]
    ndelete = int(sys.argv[2])
    output = sys.argv[3]
    
    tl = dendropy.TreeList()
    tl.read_from_path(inp, 'newick')
    
    for t in tl:
        delete_taxa(t, ndelete)
    tl.write_to_path(output, 'newick')
