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
  Clade(TaxonSet& tx, string& str);
  Clade(TaxonSet& tx, bitset<128> taxa);

  Clade(TaxonSet& tx);
  Clade(const Clade& other) :
    taxa_list(other.taxa_list), taxa(other.taxa), tx(other.tx)
  {}
  Clade& operator=(const Clade& other) {
    taxa_list = other.taxa_list;
    taxa = other.taxa;
    tx = other.tx;
    return *this;
  }
  
  string str() const;
  string newick_str(TripartitionScorer& scorer, vector<Clade>& clades);
  
  bool contains(const Clade& other) const;
  bool contains(const Taxon taxon) const;

  size_t size() { return taxa_list.size(); }
  
  static void test();
  
  void add(const Taxon taxon);
  Clade complement() const;
  Clade minus(const Clade& other) const;

  double score(TripartitionScorer& scorer, vector<Clade>& clades, unordered_set<bitset<128> >& cladetaxa);
  
  unordered_set<Taxon> taxa_list;
  bitset<128> taxa;
private:
  TaxonSet& tx;
};

struct Tripartition {
  Clade a1, a2, rest;
  Tripartition(TaxonSet& tx, Clade& clade, Clade& subclade) :
    a1(tx), a2(tx), rest(tx)
  {
    assert(clade.contains(subclade));
    a1 = clade.minus(subclade);
    a2 = subclade;
    rest = clade.complement();
  }
};


#endif // CLADE_HPP__
