#ifndef CLADE_HPP__
#define CLADE_HPP__

#include <string>
#include <vector>
#include "TaxonSet.hpp"

using namespace std;

class Clade {
public:
  Clade(TaxonSet& tx, string& str);

  string str();
  bool contains(const Clade& other);
  bool contains(const int taxon);

  static void test();
  
private:
  vector<int> taxa;
  TaxonSet& tx;
};

#endif // CLADE_HPP__
