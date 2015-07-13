#include "CladeSelector.hpp"
#include <algorithm>

void CladeSelector::run() {
  sort(clades.begin(), clades.end(), [](Clade& a, Clade& b){ return a.size() < b.size(); });
  
  for (Clade& clade : clades){
    clade.score(scorer, clades, cladetaxa);
  }
  cout << clades.back().str() << endl;
  cout << clades.back().newick_str(scorer, clades) << endl;
  cout << clades.back().score(scorer, clades, cladetaxa);
}
