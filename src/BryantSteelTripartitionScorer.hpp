#ifndef BRYANTSTEEL_TRIPARTITION_SCORER_HPP__
#define BRYANTSTEEL_TRIPARTITION_SCORER_HPP__

#include "TripartitionScorer.hpp"
#include "Logger.hpp"


class BryantSteelTripartitionScorer : public TripartitionScorer{
public:
  DEC_SCORER(BryantSteelTripartitionScorer);
  BryantSteelTripartitionScorer(TaxonSet& ts);
  virtual double score(const Tripartition<ScorableClade>& t);
private:
  unordered_map<clade_bitset, map<pair<Taxon, Taxon>, double> >  W;
  QuartetDict qd;
};



#endif
