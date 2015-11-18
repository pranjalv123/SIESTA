#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


#include <unordered_map>
#include "Clade.hpp"
#include "Quartet.hpp"

using namespace std;


namespace std {
  template <typename T, typename U> struct hash<pair<T,U> > {
    size_t operator()(pair<T, U> x) {
      return hash<T>()(x.first) ^ hash<U>()(x.second);
    }
  };
};
  
class TripartitionScorer {
public:
  virtual double score(const Tripartition& t)=0;
  double get_score(clade_bitset& clade);
  void set_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2);
  pair<clade_bitset, clade_bitset>& get_subclades(clade_bitset& clade, vector<Clade>& clades);
  TripartitionScorer(TaxonSet& ts) : ts(ts) {
    Clade ec(ts);
    score_map[ec.get_taxa()] = 0;
    subclade_map.emplace(ec.get_taxa(), make_pair(ec.get_taxa(), ec.get_taxa()));
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

class BryantSteelTripartitionScorer : public TripartitionScorer{
public:
  BryantSteelTripartitionScorer(TaxonSet& ts, QuartetDict& qd, vector<Clade>& clades);
  virtual double score(const Tripartition& t);
private:
  unordered_map<clade_bitset, map<pair<Taxon, Taxon>, double> >  W;
  QuartetDict& qd;
};


class RFTripartitionScorer : public TripartitionScorer {
public:
  RFTripartitionScorer(TaxonSet& ts, string treesfile);
  RFTripartitionScorer(TaxonSet& ts, vector<string> trees);
  virtual double score (const Tripartition& t);
  bool matches(const Tripartition& t, const Bipartition& bp);
private:
  unordered_map<Bipartition, double > clade_weights;
};

#endif
