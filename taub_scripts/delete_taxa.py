import dendropy
import numpy as np

def delete_taxa(t, seqsubset, n):
    taxa = [i.taxon for i in np.random.choice(t.leaf_nodes(), size=n, replace=False)]
    labels = set([str(i.label) for i in taxa])
    t.prune_taxa(taxa)
    modseq = []
    for i in seqsubset:
        if i.split()[0] in labels:
            modseq.append(i.replace('A', '-').replace('C', '-').replace('T', '-').replace('G', '-'))
        else:
            modseq.append(i)
#    seqsubset =[i for i in seqsubset if i.split()[0] not in labels]
    return modseq
    
    
def isboundary(s):
    splits = s.split()
    if len(splits) > 1:
        return splits[1].isdigit()
    return False

import sys
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "delete_taxa.py input seq_input ndelete..."
        exit()
    inp = sys.argv[1]
    seq_inp = sys.argv[2]
    ndelete = []
    for i in sys.argv[3:]:
        ndelete.append(int(i))
    
    ndelete.sort()
    print inp
    print seq_inp
    print ndelete
    tl = dendropy.TreeList()
    tl.read_from_path(inp, 'newick')
    startree = dendropy.treesim.star_tree(tl.taxon_set)
    seqs = open(seq_inp).readlines()

    seq_gene_starts = [i[0] for i in enumerate(seqs) if isboundary(i[1])]
    seq_gene_starts.append(-1)

    pn = 0

    for n in ndelete:
        outseq = open(seq_inp + "m" + "%02d" % n, 'w')
        out = open(inp + "m" + "%02d" % n, 'w')
        i = 0
        for i in range(len(tl)):
#        for t, s1, s2 in zip(tl, seq_gene_starts[:-1], seq_gene_starts[1:]):
            s1 = seq_gene_starts[i]
            s2 = seq_gene_starts[i+1]
            modseq = delete_taxa(tl[i], seqs[s1:s2], n - pn)
            outseq.write(''.join(modseq))

        pn = n
        tl.write_to_stream(out, 'newick')
        startree.write_to_stream(out, 'newick')
        outseq.close()
        out.close()
        
