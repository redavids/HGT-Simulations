import dendropy
import gmpy
from collections import Counter

def real_leaves(T):
    taxa = [i.taxon for i in T.leaf_nodes() if 'dummy' not in i.taxon.label]
    bitmask = 0
    for t in taxa:
        bitmask |= T.taxon_set.taxon_bitmask(t)
    return bitmask

def difference(T, t):
    diff = 0
    for e in t.get_edge_set():
        if len(edges_by_split(T, t, e.split_bitmask)) == 0:
            diff+=1
    return diff


def difference_list(T, t):
    diff = []
    for e in t.get_edge_set():
        if len(edges_by_split(T, t, e.split_bitmask)) == 0 or len(edges_by_split(T, t, t.seed_node.edge.split_bitmask ^ e.split_bitmask)) == 0:
            diff.append(min(e.split_bitmask & t.seed_node.edge.split_bitmask, (t.seed_node.edge.split_bitmask ^ e.split_bitmask) & t.seed_node.edge.split_bitmask))

    print diff

    return set(diff)

def edges_by_split(T, t, split):
    
    taxonset = t.seed_node.edge.split_bitmask
    Taxonset = T.seed_node.edge.split_bitmask


    return [e for e in T.get_edge_set() 
            if (e.split_bitmask & taxonset) == split 
#            or (e.split_bitmask & taxonset) == (taxonset ^ split) & taxonset]
            or ((Taxonset ^ e.split_bitmask) & taxonset) == split]



def pick_good_edges(T, t):

    T_splits = T.split_edges
    t_splits = t.split_edges

    all_t_bitmask = t.seed_node.edge.split_bitmask
    all_T_bitmask = T.seed_node.edge.split_bitmask

    TA_splits = {}
    for split in T_splits:
        TA_splits[split & all_t_bitmask] = T_splits[split]
        TA_splits[t.taxon_set.complement_split_bitmask(split) & all_t_bitmask] = T_splits[split]
    
    matching_splits = []
    for split in t_splits:
        if split in TA_splits:
            if split != all_t_bitmask and split != 0:
                matching_splits.append(split)
    if all_t_bitmask in TA_splits:
        matching_splits.append(all_t_bitmask)

    return matching_splits

def get_subtrees(T, t, edges, s1_bitmask, s2_bitmask, tx1, tx2):
    T.taxon_set.add(tx1)
    T.taxon_set.add(tx2)

    T.reroot_at_edge(edges[0], update_splits = True)

#    try:
    n1 = T.mrca(split_bitmask = s1_bitmask)
    S1_bitmask = n1.edge.split_bitmask
    n2 = T.mrca(split_bitmask = s2_bitmask)
    S2_bitmask = n2.edge.split_bitmask

#    except ValueError:
#        T.print_plot()
#        print s1_bitmask
#        print s2_bitmask
#        print sorted([i.split_bitmask for i in T.get_edge_set()])[::-1]
#        pass

    

    T1 = dendropy.Tree(T)
    T2 = dendropy.Tree(T)
    Path = dendropy.Tree(T)

    T1.retain_taxa(T1.taxon_set.split_taxa_list(S1_bitmask))
    T2.retain_taxa(T2.taxon_set.split_taxa_list(S2_bitmask))

    if T1.seed_node.taxon != None:
        T1.seed_node.add_child(dendropy.Node(taxon=T1.seed_node.taxon))
        T1.seed_node.taxon = None
    n1 = dendropy.Node(taxon=tx1)
    T1.seed_node.add_child(n1)

    if T2.seed_node.taxon != None:
        T2.seed_node.add_child(dendropy.Node(taxon=T1.seed_node.taxon))
        T2.seed_node.taxon = None
    n2 = dendropy.Node(taxon=tx2)
    T2.seed_node.add_child(n2)
    
    
    pn1 = Path.mrca(split_bitmask = s1_bitmask)
    pn2 = Path.mrca(split_bitmask = s2_bitmask)


    pn1.add_child(dendropy.Node(taxon=tx1))
    Path.update_splits()
    pn2.add_child(dendropy.Node(taxon=tx2))
    Path.update_splits()


    Path.prune_taxa(Path.taxon_set.split_taxa_list(S1_bitmask))
    Path.update_splits()
    Path.prune_taxa(Path.taxon_set.split_taxa_list(S2_bitmask))
    Path.update_splits()


    return T1, T2, Path, tx1, tx2, pn1, pn2, n1, n2


def root_t_correctly(T, t, s1, s2):
    print "used special case"
    initial_nleaf = len(T.leaf_nodes())
    T_orig = dendropy.Tree(T)
    t_orig = dendropy.Tree(t)
    s2 = (T.seed_node.edge.split_bitmask ^ s1) & T.seed_node.edge.split_bitmask
#    print [i for i in T.split_edges if i == s1]
    edges = edges_by_split(T, t, s1)
#    print edges
    edge = edges[0]
#    T.reroot_at_edge(edge)

    tx1 = dendropy.Taxon(label='dummy-3')
    T.taxon_set.add(tx1)

    n2 = None

    if T.mrca(split_bitmask=s2) == T.seed_node:
        n2 = T.mrca(split_bitmask = s1)
        n2.add_child(dendropy.Node(taxon=tx1))
        T.update_splits()

#        T.reseed_at(T.mrca(split_bitmask=s1))
#        T.update_splits()
    elif T.mrca(split_bitmask=s1) == T.seed_node:
        n2 = T.mrca(split_bitmask = s2)
        n2.add_child(dendropy.Node(taxon=tx1))
        T.update_splits()

#        T.reseed_at(T.mrca(split_bitmask=s2))
#        T.update_splits()
    else:
        print "huh?"

    
    tmask = n2.child_nodes()[0].edge.split_bitmask
    TA = dendropy.Tree(T)
    TA.retain_taxa(T.taxon_set.split_taxa_list(s1))
    T.prune_taxa(T.taxon_set.split_taxa_list(s1))
    T.update_splits()
    TA.update_splits()
    t.update_splits()

    split = TA.seed_node.child_nodes()[0].edge.split_bitmask & t.seed_node.edge.split_bitmask
    
    print TA.taxon_set.split_taxa_list(split)
    TA.print_plot()
    if t.mrca(split_bitmask = split ^ t.seed_node.edge.split_bitmask) == t.mrca(split_bitmask = split):
        t.reroot_at_edge(t.seed_node.child_nodes()[0].edge)
    elif t.mrca(split_bitmask=split) in [None, t.seed_node]:
        t.reroot_at_node(t.mrca(split_bitmask = split ^ t.seed_node.edge.split_bitmask))
    else: 
        n = t.mrca(split_bitmask=split)
        if n.is_leaf:
            t.reroot_at_edge(n.incident_edges()[0])
        else:
            t.reroot_at_node(n)
        
    

    t.update_splits()
    T.update_splits()
    
    if len(T.leaf_nodes()) == 1:
        print sorted([i.taxon.label for i in TA.leaf_nodes()])
        print sorted([i.taxon.label for i in t.leaf_nodes()])
        print sorted([i.taxon.label for i in T.leaf_nodes()])
        T.clone_from(t)
        return
    
#    T.print_plot(plot_metric='level', show_internal_node_ids=True)
#    t.print_plot(plot_metric='level', show_internal_node_ids=True)
    n1 = [i for i in T.nodes() if i.taxon == tx1][0]
    if n1.parent_node == None:
        T.print_plot(plot_metric='level')
        TA.print_plot(plot_metric='level')
        t.print_plot(plot_metric='level')
        T.reroot_at_edge(n1.incident_edges()[0])
    n1.parent_node.add_child(t.seed_node)

        

    T.update_splits()
    t.update_splits()
    T.prune_taxa([tx1])
    if False:
 #   if not (len(T.leaf_nodes()) == initial_nleaf):
        print [i.label for i in T.taxon_set.split_taxa_list(s1)]
        print len(T.leaf_nodes()), initial_nleaf
        print n1
        T.print_plot(plot_metric='level', show_internal_node_ids=True)
        T_orig.print_plot()
        t.print_plot()
        t_orig.print_plot(plot_metric='level')
#        assert(False)
    
def match_subtrees_(T, t):    

    T.update_splits()
    t.update_splits()

    t.deroot()
    T.deroot()
    
    goodedges = pick_good_edges(T, t) 

    edgeswithsubtrees = [e for e in goodedges if len(edges_by_split(T, t, e)) > 1]

    if real_leaves(T) == real_leaves(t):
        return t

    if len(edgeswithsubtrees) == 0:
        return t

    split = edgeswithsubtrees[0]


    s1 = split
    s2 = t.seed_node.edge.split_bitmask ^ s1

    if s2 == 0:
        s2 = T.seed_node.edge.split_bitmask ^ s1
        root_t_correctly(T, t, s1, s2)
        return T


    tx1 = dendropy.Taxon(label='dummy-1')
    tx2 = dendropy.Taxon(label='dummy-2')  

    T_edges = edges_by_split(T, t, split)
    T1, T2, Path, Tx1, Tx2, PN1, PN2, N1, N2 = get_subtrees(T, t, T_edges, s1, s2, tx1, tx2)        
    
    t_edge = edges_by_split(t, t, split)
    t1, t2, path, tx1, tx2, pn1, pn2, n1, n2 = get_subtrees(t, t, t_edge, s1, s2, tx1, tx2)

    t1 = match_subtrees_(T1, t1)
    t2 = match_subtrees_(T2, t2)

    n1 = [i for i in t1.nodes() if i.taxon == tx1][0]
    n2 = [i for i in t2.nodes() if i.taxon == tx2][0]

    PN1 = [i for i in Path.nodes() if i.taxon == tx1][0]
    PN2 = [i for i in Path.nodes() if i.taxon == tx2][0]

    t1.reseed_at(n1.parent_node)
    t2.reseed_at(n2.parent_node)

    PN1.add_child(n1.parent_node)
    PN2.add_child(n2.parent_node)
    
    return Path
        
def match_subtrees(T, t):
    Tx = dendropy.Tree()
    tx = dendropy.Tree()
    Tx.clone_from(T)
    tx.clone_from(t)
    tr = match_subtrees_(Tx, tx)
    
    tr.update_splits()
    
    labels = [i for i in tr.taxon_set if 'dummy' in i.label]


    tr.prune_taxa(labels)

    return tr
                


import sys
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "python supertree.py tree_file output_file"
    tree_file = sys.argv[1]
    output_file = open(sys.argv[2], 'w')
    big_trees = []
    little_trees = []
    completed_trees = []
    all_tree_list = dendropy.TreeList()
    
    all_tree_list.read_from_path(tree_file, 'newick')

    for tree in all_tree_list:
        print len(tree.leaf_nodes())
        if len(tree.leaf_nodes()) == len(tree.taxon_set):
            big_trees.append(tree)
        else:
            little_trees.append(tree)
    k = 0
    ta = dendropy.Tree()
    for T in big_trees:
        for t in little_trees:
            k += 1
#            if k != 14:
#                continue
            tx = match_subtrees(T, t)
            output_file.write(str(tx) + ';' + '\n')
            completed_trees.append(tx)
            ta.clone_from(T)
            ta.update_splits()
            T.update_splits()
            tx.update_splits()
            t.update_splits()
            
#            tr = dendropy.Tree(tx)
#            tr.update_splits()
#            tr.retain_taxa(t.taxon_set.split_taxa_list(t.seed_node.edge.split_bitmask))

            T.deroot()
            print "=="
            print difference(tx, t), '= 0'
            print difference(tx, T), '=', difference(T, t)
            print len(T.leaf_nodes()), '=', len(tx.leaf_nodes())
            print len(T.get_edge_set()), '=', len(tx.get_edge_set())
            dl1 =  difference_list(T, t)
            dl2 = difference_list(tx, T)

    output_file.close()
    
