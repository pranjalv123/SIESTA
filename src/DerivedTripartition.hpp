#ifndef DERIVED_TRIPARTITION_SCORER_HPP__
#define DERIVED_TRIPARTITION_SCORER_HPP__


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

class BryantSteelTripartitionScorer : public TripartitionScorer{
public:
  DEC_SCORER(BryantSteelTripartitionScorer);
  BryantSteelTripartitionScorer(TaxonSet& ts);
  virtual double score(const Tripartition& t);
private:
  unordered_map<clade_bitset, map<pair<Taxon, Taxon>, double> >  W;
  QuartetDict qd;
};


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

#endif // DERIVED_TRIPARTITION_SCORER_HPP__
