#ifndef CLADE_SELECTOR_HPP__
#define CLADE_SELECTOR_HPP__

class CladeTreeNode {
  
};

class CladeSelector {
  TripartitionScorer& scorer;
  vector<Clade>& clades;
  TaxonSet& ts;
  
public:
  CladeSelector() {}

  void run();
};

#endif
