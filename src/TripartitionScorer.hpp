#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


#include <unordered_map>
#include "Clade.hpp"
#include "Quartet.hpp"

class TripartitionScorer {
public:
  virtual double score(const Tripartition& t)=0;
  double get_score(bitset<128>& clade);
  void set_score(bitset<128>& clade, double score, Clade& a1, Clade& a2);
  pair<Clade, Clade>& get_subclades(bitset<128>& clade, vector<Clade>& clades);
  TripartitionScorer(TaxonSet& ts) : ts(ts) {}
 
private:
  TaxonSet& ts;
  unordered_map <bitset<128>, double> score_map;
  unordered_map <bitset<128>, pair<Clade, Clade> > subclade_map;
};

class DPTripartitionScorer : public TripartitionScorer{
public:
  DPTripartitionScorer(TaxonSet& ts, QuartetDict& qd) : TripartitionScorer(ts), qd(qd) {}
  virtual double score(const Tripartition& t);
private:
  QuartetDict& qd;
};

#endif
