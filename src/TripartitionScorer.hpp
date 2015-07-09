#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


struct Tripartition {
  Clade a1, a2, rest;
  Tripartition(TaxonSet& tx, Clade& clade, Clade& subclade) {
    a1 = clade.minus(subclade);
    a2 = subclade;
    rest = clade.complement();
  }
};

class TripartitionScorer {
public:
  virtual double score(const Tripartition& t);
};

class DPTripartitionScorer : TripartitionScorer{
public:
  DPTripartitionScorer(TaxonSet& ts, QuartetDict& qd) : ts(ts), qd(qd) {}
  virtual double score(const Tripartition& t);
private:
  TaxonSet& ts;
  QuartetDict& qd;
};

#endif
