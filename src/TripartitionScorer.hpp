#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


#include <unordered_map>
#include <map>
#include <Logger.hpp>
#include <Clade.hpp>
#include <Quartet.hpp>

using namespace std;


namespace std {
  template <typename T, typename U> struct hash<pair<T,U> > {
    size_t operator()(pair<T, U> x) {
      return hash<T>()(x.first) ^ hash<U>()(x.second);
    }
  };
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
