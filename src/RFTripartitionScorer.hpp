#ifndef RF_TRIPARTITION_SCORER_HPP__
#define RF_TRIPARTITION_SCORER_HPP__

#include "TripartitionScorer.hpp"

class RFTripartitionScorer : public TripartitionScorer {
public:
  DEC_SCORER(RFTripartitionScorer);
  RFTripartitionScorer(TaxonSet& ts);
  int addSourceTree(string tree);
  virtual double score (const Tripartition& t);
  bool matches(const Tripartition& t, const Bipartition& bp);
  int total_weight;
  int n_trees;
  virtual double adjust_final_score(double score);
private:
  unordered_map<Bipartition, double > clade_weights;
};



#endif // RF_TRIPARTITION_SCORER_HPP__
