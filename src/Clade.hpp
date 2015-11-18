#ifndef CLADE_HPP__
#define CLADE_HPP__

#include <string>
#include <vector>
#include <iostream>
#include <bitset>
#include <unordered_set>
#include <cassert>
#include <unordered_set>
#include <string.h>
#include "TaxonSet.hpp"


class TripartitionScorer;
using namespace std;

class Clade {
private:


public:
  clade_bitset taxa;
  TaxonSet& ts;
  
  Clade(TaxonSet& ts, string& str);
  Clade(TaxonSet& ts, clade_bitset& taxa);
  Clade(TaxonSet& ts, unordered_set<Taxon>& taxa);
  Clade(TaxonSet& ts);
  Clade(const Clade& other);
  
  
  Clade& operator=(const Clade& other);
  bool operator==(const Clade& other) const;
  
  
  string str() const;
  string newick_str(TripartitionScorer& scorer, vector<Clade>& clades);
  
  bool contains(const Clade& other) const;
  bool contains(const Taxon taxon) const;

  static void test();
  
  void add(const Taxon taxon);
  Clade complement() const;
  Clade minus(const Clade& other) const;

  int size() const;
  const clade_bitset& get_taxa() const {return taxa;}
  
  double score(TripartitionScorer& scorer, vector<Clade>& clades, unordered_set<clade_bitset>& cladetaxa);

  BVFIterator begin() const {
    return taxa.begin();
  }
  
  BVFIterator end() const {
    return taxa.end();
  }

  void do_swap(Clade& other);
  size_t hash() const { return taxa.hash(); }
};

struct Tripartition {
  Clade a1, a2, rest;
  Tripartition(TaxonSet& ts, Clade& clade, Clade& subclade);
};

struct Bipartition {
  Clade a1, a2;
  Bipartition(const Clade& clade1, const Clade& clade2) :
    a1(clade1),
    a2(clade2)
  {}
  size_t hash() const { return a1.taxa.hash() ^ a2.taxa.hash(); }
  bool operator==(const Bipartition& other) const {
    return ((a1 == other.a1) && (a2 == other.a2)) || ((a2 == other.a1) && (a1 == other.a2));
  }
};


namespace std {
  template <> struct hash<Clade> {
    size_t operator()(const Clade& bvf) const {
      return bvf.hash();
    }
  };
    template <> struct hash<Bipartition> {
    size_t operator()(const Bipartition& bp) const {
      return bp.hash();
    }
  };
}


#endif // CLADE_HPP__
