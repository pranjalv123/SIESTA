#ifndef CLADE_HPP__
#define CLADE_HPP__

#include <string>
#include <vector>
#include <bitset>
#include <unordered_set>
#include <cassert>
#include <unordered_set>
#include "TaxonSet.hpp"


class TripartitionScorer;
using namespace std;

class Clade {
public:
  Clade(TaxonSet& ts, string& str);
  Clade(TaxonSet& ts, clade_bitset taxa);
  Clade(TaxonSet& ts);
  Clade(const Clade& other);
  
  Clade& operator=(const Clade& other);
  
  string str() const;
  string newick_str(TripartitionScorer& scorer, vector<Clade>& clades);
  
  bool contains(const Clade& other) const;
  bool contains(const Taxon taxon) const;

  size_t size();
  
  static void test();
  
  void add(const Taxon taxon);
  Clade complement() const;
  Clade minus(const Clade& other) const;

  double score(TripartitionScorer& scorer, vector<Clade>& clades, unordered_set<clade_bitset>& cladetaxa);
  
  unordered_set<Taxon> taxa_list;
  clade_bitset taxa;

  TaxonSet& ts;
};

struct Tripartition {
  Clade a1, a2, rest;
  Tripartition(TaxonSet& ts, Clade& clade, Clade& subclade);
};


#endif // CLADE_HPP__
