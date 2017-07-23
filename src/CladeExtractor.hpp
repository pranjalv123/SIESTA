#ifndef CLADE_EXTRACTOR__
#define CLADE_EXTRACTOR__
#include <vector>
#include <string>
#include <Clade.hpp>
#include <util/Options.hpp>

using namespace std;

class CladeExtractor {
public:
  static unordered_set<Clade> extract(TaxonSet& ts, const string& tree, unordered_set<Taxon>& tree_taxa);
  static void test();

  static TaxonSet& get_taxonset();
  static unordered_set<Clade>& get_clades();
  static unordered_set<clade_bitset>& get_cladetaxa();
private:
  static void get_from_cl();
  static unordered_set<Clade> cl_clades;
  static unordered_set<clade_bitset > cl_cladetaxa;
  static TaxonSet* ts;
};
#endif
