#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


#include <unordered_map>
#include <map>
#include <sstream>
#include <Clade.hpp>
#include <util/Logger.hpp>
#include <Quartet.hpp>
#include "ScorableClade.hpp"

class ScorableClade;
class ScorableTripartition;

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
  //Returns the score of a tripartition
  virtual double score(const Tripartition<ScorableClade>& t)=0;
  virtual double adjust_final_score(double score);

  // these are used internally by the DP algorithm
  double get_score(clade_bitset& clade);
  void set_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2);
  void add_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2);
  void clear_scores(clade_bitset& clade);
  pair<clade_bitset, clade_bitset>& get_subclades(clade_bitset& clade, vector<ScorableClade>& clades);
  vector<pair<clade_bitset, clade_bitset> >& get_subclade_lists(clade_bitset& clade);
    
  TripartitionScorer(TaxonSet& ts) : ts(ts) {
    Clade ec(ts);
    score_map[ec.get_taxa()] = 0;
    subclade_map.emplace(ec.get_taxa(), make_pair(ec.get_taxa(), ec.get_taxa()));
  }
protected:
  TaxonSet& ts;
private:
  unordered_map <clade_bitset, double> score_map;
  unordered_map <clade_bitset, pair<clade_bitset, clade_bitset> > subclade_map;
  unordered_map <clade_bitset, vector<pair<clade_bitset, clade_bitset> > > subclade_list_map;

};




//below is fun internal stuff 


template <typename T> TripartitionScorer* createT(TaxonSet& ts) { return new T(ts); }

struct TripartitionScorerFactory {
  typedef map<string, TripartitionScorer*(*)(TaxonSet& ts)> map_type;

  static TripartitionScorer * createInstance(string const& s, TaxonSet& ts) {
    map_type::iterator it = getMap()->find(s);
    if(it == getMap()->end()) {
      map_type* mp = getMap();
      INFO << "Valid criteria:" << endl;
      for (auto& i : *mp) {
	INFO << i.first << endl;
      }
      return 0;
    }
    return it->second(ts);      
  }
  static string instanceList() {
    stringstream ss;
    for (auto const& it : *getMap()) {
      ss << it.first << "\n";
    }
    return ss.str();
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


#endif
