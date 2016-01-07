#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


#include <unordered_map>
#include <map>
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
  virtual double adjust_final_score(double score);
protected:
  TaxonSet& ts;
private:
  unordered_map <clade_bitset, double> score_map;
  unordered_map <clade_bitset, pair<clade_bitset, clade_bitset> > subclade_map;
};


template <typename T> TripartitionScorer* createT(TaxonSet& ts) { return new T(ts); }

struct TripartitionScorerFactory {
  typedef map<string, TripartitionScorer*(*)(TaxonSet& ts)> map_type;

  static TripartitionScorer * createInstance(string const& s, TaxonSet& ts) {
    map_type::iterator it = getMap()->find(s);
    if(it == getMap()->end())
      return 0;
    return it->second(ts);      
  }
 protected:
  static map_type* getMap() {
    if (!mymap) { mymap = new map_type; }
    return mymap;
  }
 private:
  static map_type *mymap;
};


template<typename T>
struct DerivedTPScorer : TripartitionScorerFactory {
  DerivedTPScorer(string const& s) {
    getMap() -> insert(make_pair(s, &createT<T>));
  }
};

#define DEC_SCORER(NAME)			\
    static DerivedTPScorer<NAME> reg

#define DEF_SCORER(NAME) \
    DerivedTPScorer<NAME> NAME::reg(#NAME)


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


#endif
