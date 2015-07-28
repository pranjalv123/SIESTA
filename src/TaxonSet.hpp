#ifndef TAXONSET_HPP__
#define TAXONSET_HPP__

#include <string>
#include <vector>
#include <bitset>
#include <map>
#include <sstream>
using namespace std;

typedef int Taxon;
typedef bitset<40> clade_bitset;

class TaxonSet {
private:
  vector<string> taxa;
  map<string, Taxon> index;
public:
  clade_bitset taxa_bs;
  Taxon operator[](const string& str) {
    return add(str);
  }
  const string& operator[](const Taxon i) const {
    return taxa.at(i);
  }
  int size() const {
    return taxa.size();
  }
  Taxon add(const string& str) {
    if (index.count(str)) {
      return index[str];
    }
    int i = taxa.size();
    taxa.push_back(str);
    index[str] = i;
    taxa_bs.set(i);
    return i;
  }  
  string str() {
    stringstream ss;
    for (int i = 0; i < taxa.size(); i++) {
      ss << i << "\t" << taxa[i] << endl;
    }
    return ss.str();
  }
};

#endif // TAXONSET_HPP__
