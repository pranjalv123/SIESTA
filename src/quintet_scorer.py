import dendropy
import argparse
from collections import defaultdict

taxon_namespace = None
quintet_weights = defaultdict(lambda: 0)


class Quintet:
    def __init__(self, t):
        t.update_bipartitions()
        self.a = None
        self.b = None
        self.c = None
        self.taxon_namespace = t.taxon_namespace
        
        childs = t.seed_node.child_nodes()
        if len(childs[0].leaf_nodes()) < 3:
            self.c = childs[0].leafset_bitmask
            self.a = childs[1].child_nodes()[0].leafset_bitmask
            self.b = childs[1].child_nodes()[1].leafset_bitmask
        if len(childs[1].leaf_nodes()) < 3:
            self.c = childs[1].leafset_bitmask
            self.a = childs[0].child_nodes()[0].leafset_bitmask
            self.b = childs[0].child_nodes()[1].leafset_bitmask


    def __repr__(self):
        return '\n'.join([str(self.taxon_namespace.bitmask_taxa_list(self.a)),str(self.taxon_namespace.bitmask_taxa_list(self.b)), str(self.taxon_namespace.bitmask_taxa_list(self.c))]) + '\n'

    def matches(self, a1, a2, rest):
        if ((a1 & self.a) == self.a) and ((a2 & self.b) == self.b) and ((rest & self.c) == self.c):
            return True
        if ((a2 & self.a) == self.a) and ((a1 & self.b) == self.b) and ((rest & self.c) == self.c):
            return True
        
        return False
        
import sys

def init(ts, opts):
    global taxon_namespace, quintet_weights

    parser = argparse.ArgumentParser(prog="wASTRAL")
    parser.add_argument('-q', '--quintets')
    args, _ = parser.parse_known_args(opts)

    taxon_namespace = ts
    
    for line in open(args.quintets).readlines():
        quintet, weight = line.split(':')
        weight = float(weight)
        quintet = dendropy.Tree.get_from_string(quintet + ';', 'newick', taxon_namespace=taxon_namespace, preserve_underscores=True, rooting='force-rooted')
        quintet.update_bipartitions()
        
            
        quintet_weights[Quintet(quintet)] = weight
    print quintet_weights
    
def score(a1, a2, rest):
    global taxon_namespace, quintet_weights
    weight = 0.0

    for q,w in quintet_weights.items():
        if q.matches(a1, a2, rest):
            weight += w

#    print taxon_namespace.bitmask_taxa_list(a1), taxon_namespace.bitmask_taxa_list(a2), taxon_namespace.bitmask_taxa_list(rest), weight
    return weight

def adjust(n):
    return n
