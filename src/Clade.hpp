#ifndef CLADE_HPP__
#define CLADE_HPP__

#include <string>
#include <vector>
#include <bitset>
#include <unordered_set>
#include "TaxonSet.hpp"

using namespace std;

class Clade {
public:
  Clade(TaxonSet& tx, string& str);
  Clade(TaxonSet& tx);
  
  
  string str();
  bool contains(const Clade& other) const;
  bool contains(const int taxon) const;
  
  static void test();
  
  void add(const Taxon taxon);
  Clade complement() const;
  Clade minus(const Clade& other) const;

  double score(TripartitionScorer& scorer);
  
  double value;
  Clade* best_subclade;
  unordered_set<int> taxa_list;
  bitset<128> taxa;
private:
  TaxonSet& tx;
};

#endif // CLADE_HPP__
