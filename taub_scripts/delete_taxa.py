import dendropy
import numpy as np

def delete_taxa(t, n):
    t.prune_taxa([i.taxon for i in np.random.choice(t.leaf_nodes(), size=n)])
    

import sys
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "delete_taxa.py input ndelete..."
        exit()
    inp = sys.argv[1]
    ndelete = []
    for i in sys.argv[2:]:
        ndelete.append(int(i))
    
    ndelete.sort()

    tl = dendropy.TreeList()
    tl.read_from_path(inp, 'newick')
    tl.append(dendropy.treesim.star_tree(tl.taxon_set))
    pn = 0

    for n in ndelete:
        for t in tl:
            if len(t.get_edge_set()) > len(t.leaf_nodes()) + 2:
                delete_taxa(t, n - pn)
        pn = n
        tl.write_to_path(inp + "m" + "%02d" % n, 'newick')
