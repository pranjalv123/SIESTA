import dendropy
import argparse
from collections import defaultdict

taxon_namespace = None
clade_weights = defaultdict(lambda: 0)
genetrees = None

import sys

def init(ts, opts):
    global genetrees, taxon_namespace, clade_weights, ts_size
    
    taxon_namespace = ts

    ts_size = len(ts)
    
    parser = argparse.ArgumentParser(prog="wASTRAL")

    parser.add_argument('-g', '--genetrees')
    
    args, _ = parser.parse_known_args(opts)

    genetrees = dendropy.TreeList.get_from_path(args.genetrees, 'newick', taxon_namespace = ts, preserve_underscores=True)

    for i in enumerate(ts):
        print i
    
    for t in genetrees:
        t.update_bipartitions()
        for b in t.encode_bipartitions():
            x = (b._split_bitmask & b._tree_leafset_bitmask, (b._tree_leafset_bitmask ^ b._split_bitmask)& b._tree_leafset_bitmask)
            if min(x) and max(x) and not b.is_trivial():
                clade_weights[(min(x), max(x))] += 1
    assert(len(ts) == ts_size)
def matches(a1, a2, b1, b2):
    if ((a1 & b1) == b1) and ((a1 & b2) == 0) and ((a2 & b2) != 0):
        return True
    if ((a2 & b1) == b1) and ((a2 & b2) == 0) and ((a1 & b2) != 0):
        return True
    if ((a1 & b2) == b2) and ((a1 & b1) == 0) and ((a2 & b1) != 0):
        return True
    if ((a2 & b2) == b2) and ((a2 & b1) == 0) and ((a1 & b1) != 0):
        return True
    return False
    
def score(a1, a2, rest):
    assert(len(taxon_namespace) == ts_size)
    weight = 0.0

    for k in clade_weights:
        if matches(a1, a2, k[0], k[1]):
            weight += clade_weights[k]
    return weight

def adjust(n):
    
    print "Total weight", sum (clade_weights.values())
    print "Raw score", n
    return (sum(clade_weights.values()) - n) * 2

