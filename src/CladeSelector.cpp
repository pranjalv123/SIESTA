#include "CladeSelector.hpp"
#include <algorithm>

CladeSelector::run() {
  sort(clades.begin(), clades.end(), [](Clade& a, Clade& b){ return a.size() < b.size(); });
  
  for (Clade& clade : clades) {
    if (clade.size() == 1) {
      clade.value = 0;
    }

    if (clade.size() == 2) {
      Clade c1;
      Clade c2;
      c1.add(*clade.taxa_list.begin());
      c2.add(*clade.taxa_list.end());
      clade.value = scorer.score(Tripartition(a.tx, c1, c2)); 
    }

    if (clade.size() > 2) {
      
    }
  }
  
}
