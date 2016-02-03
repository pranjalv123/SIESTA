#ifndef DP_TRIPARTITION_SCORER_HPP__
#define DP_TRIPARTITION_SCORER_HPP__

#include "TripartitionScorer.hpp"



//this is the kind of thing you need to do to implement a tripartition
//scorer - make a constructor and a score function, you also need to
//call DEC_SCORER(name of scorer) in the body of the class, and
//DEF_SCORER(name of scorer) in the .cpp file.

class DPTripartitionScorer : public TripartitionScorer{
public:
  DEC_SCORER(DPTripartitionScorer);
  DPTripartitionScorer(TaxonSet& ts);
  virtual double score(const Tripartition& t);
private:
  QuartetDict& qd;
};



#endif
