#ifndef CLADE_SELECTOR_HPP__
#define CLADE_SELECTOR_HPP__

#include "TripartitionScorer.hpp"
#include <unordered_set>

class CladeTreeNode {
  
};

class CladeSelector {
  TripartitionScorer& scorer;
  vector<Clade>& clades;
  unordered_set<bitset<128> >& cladetaxa;
  TaxonSet& ts;
  
public:
  CladeSelector(TaxonSet& ts, TripartitionScorer& scorer, vector<Clade>& clades, unordered_set<bitset<128> >& cladetaxa):
    scorer(scorer),
    clades(clades),
    cladetaxa(cladetaxa),
    ts(ts)
  {}

  void run();
};

#endif
