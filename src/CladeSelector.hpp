#ifndef CLADE_SELECTOR_HPP__
#define CLADE_SELECTOR_HPP__

#include "TripartitionScorer.hpp"
#include "ScorableClade.hpp"
#include <unordered_set>

class CladeSelector {
protected:
  TripartitionScorer& scorer;
  vector<ScorableClade>& clades;
  unordered_set<clade_bitset >& cladetaxa;
  TaxonSet& ts;

  
public:
  CladeSelector(TaxonSet& ts, TripartitionScorer& scorer, vector<ScorableClade>& clades, unordered_set<clade_bitset >& cladetaxa):
    scorer(scorer),
    clades(clades),
    cladetaxa(cladetaxa),
    ts(ts)
  {}

  virtual double run(bool invert, twod_mat* mat) = 0;
  string newick_tree;
};

class BasicCladeSelector : public CladeSelector {
public:
  using CladeSelector::CladeSelector;
  virtual double run(bool invert, twod_mat* mat);
};

#endif
