#ifndef FASTRF_TRIPARTITION_SCORER_HPP__
#define FASTRF_TRIPARTITION_SCORER_HPP__

#include "Clade.hpp"
#include "BitVector.hpp"
#include "TripartitionScorer.hpp"

class FastRFTripartitionScorer : public TripartitionScorer {
public:
  DEC_SCORER(FastRFTripartitionScorer);
  FastRFTripartitionScorer(TaxonSet& ts);
  int addSourceTree(string tree);
  virtual double score (const Tripartition<ScorableClade>& t);
  bool matches(const Tripartition<ScorableClade>& t, const Bipartition& bp);
  int total_weight;
  int n_trees;
  virtual double adjust_final_score(double score);
private:
  unordered_map<Bipartition, double > clade_weights;

  vector<const Bipartition*> bipartitions;
  vector<double> weights;

  
  unordered_map<Clade, BitVectorFixed> possibleA1; //contains all of a1, none of a2
  unordered_map<Clade, BitVectorFixed> possibleB1; //contains at least one of a1, none of a2
  unordered_map<Clade, BitVectorFixed> possibleA2; //contains all of a2, none of a1
  unordered_map<Clade, BitVectorFixed> possibleB2; //contains at least one of a2, none of a1
};



#endif // RF_TRIPARTITION_SCORER_HPP__
