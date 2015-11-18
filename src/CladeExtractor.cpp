#include "CladeExtractor.hpp"
#include "TaxonSet.hpp"
#include <vector>
#include <string>
#include <cstring>
#include <cassert>

unordered_set<Clade> CladeExtractor::extract(TaxonSet& ts, const string& tree_c, unordered_set<Taxon>& tree_taxa) {
  string tree(tree_c);
  vector<size_t> active;
  vector<Clade> clades;

  char *n_cp = 0, *n_op = 0, *n_cm = 0; //next close paren, open paren, comma

  char* current = &(tree[0]);
  char* end = &(tree[tree.length()]);

  while (current < end) {
    //    cout << current << endl;
    
    n_cp = strchr(current, ')');
    n_op = strchr(current, '(');
    n_cm = strchr(current, ',');

    if (n_cp == 0 && n_op == 0 && n_cm == 0) {
      return unordered_set<Clade>(clades.begin(), clades.end());
    }
    if (n_cp == 0) 
      n_cp = end;
    if (n_op == 0)
      n_op = end;
    if (n_cm == 0)
      n_cm = end;
    
    if (n_cp < n_op && n_cp < n_cm) { //the next delim is a close paren
      if(current < n_cp) { //there is a token 
	string token(current, n_cp - current);
	Taxon taxon = ts[token];
	tree_taxa.insert(taxon);
	cout << "found " << token << endl;
	for (int a : active) {
	  clades.at(a).add(taxon);
	}
      }
      //now a clade is complete
      active.pop_back();
      current = n_cp + 1;
    } else if (n_op < n_cp && n_op < n_cm) { //the next delim is an open paren
      if(current < n_op) { //there is a token 
	string token(current, n_op - current);
	Taxon taxon = ts[token];
	tree_taxa.insert(taxon);
	//	cout << "found " << token << endl;
	for (int a : active)
	  clades.at(a).add(taxon);

      }
      //now add a new clade, and then we'll add taxa to it
      clades.emplace_back(ts);
      active.push_back(clades.size() - 1);
      current = n_op + 1;
    } else { //the next delim is a comma
      if(current < n_cm) { //there is a token
	string token(current, n_cm - current);
	Taxon taxon = ts[token];
	tree_taxa.insert(taxon);
	//	cout << "found " << token << endl;
	for (int a : active)
	  clades.at(a).add(taxon);
      }
      current = n_cm + 1;
    }
    
  }  
  return unordered_set<Clade>(clades.begin(), clades.end());
}


void CladeExtractor::test() {
  cout << "Testing CladeExtractor..." << endl;
  TaxonSet ts("");
  string tree("(((a,b),(c,d)),e,(f,g,(h,i)));");
  cout << tree << endl;
  unordered_set<Taxon> tree_taxa;
  unordered_set<Clade> clades = CladeExtractor::extract(ts, tree, tree_taxa);
  for (const Clade& c : clades) {
    cout << c.str() << endl;
  }
  cout << "Done!" << endl;
}
