#ifndef CLADE_EXTRACTOR__
#define CLADE_EXTRACTOR__
#include <vector>
#include <string>
#include "Clade.hpp"

using namespace std;

class CladeExtractor {
public:
  static unordered_set<Clade> extract(TaxonSet& ts, const string& tree, unordered_set<Taxon>& tree_taxa);
  static void test();
};

#endif
