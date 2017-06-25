#ifndef BRYANTSTEEL_TRIPARTITION_SCORER_HPP__
#define BRYANTSTEEL_TRIPARTITION_SCORER_HPP__

#include "TripartitionScorer.hpp"
#include "Logger.hpp"


class ASTRALTripartitionScorer : public TripartitionScorer{
public:
  DEC_SCORER(ASTRALTripartitionScorer);
  ASTRALTripartitionScorer(TaxonSet& ts);
  virtual double score(const Tripartition<ScorableClade>& t);
  virtual double adjust_final_score(double score);
private:
  vector<vector<Taxon> > postorders;
  vector<Clade > gttaxa;
};



#endif
