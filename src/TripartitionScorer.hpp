#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


#include <unordered_map>
#include "Clade.hpp"
#include "Quartet.hpp"

class TripartitionScorer {
public:
  virtual double score(const Tripartition& t)=0;
  double get_score(clade_bitset& clade);
  void set_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2);
  pair<clade_bitset, clade_bitset>& get_subclades(clade_bitset& clade, vector<Clade>& clades);
  TripartitionScorer(TaxonSet& ts) : ts(ts) {
    Clade ec(ts);
    score_map[ec.taxa] = 0;
    subclade_map[ec.taxa] = make_pair(ec.taxa, ec.taxa);
  }
 
private:
  TaxonSet& ts;
  unordered_map <clade_bitset, double> score_map;
  unordered_map <clade_bitset, pair<clade_bitset, clade_bitset> > subclade_map;
};

class DPTripartitionScorer : public TripartitionScorer{
public:
  DPTripartitionScorer(TaxonSet& ts, QuartetDict& qd) : TripartitionScorer(ts), qd(qd) {}
  virtual double score(const Tripartition& t);
private:
  QuartetDict& qd;
};

#endif
