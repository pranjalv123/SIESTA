#include "CladeSelector.hpp"
#include <algorithm>





void CladeSelector::run(bool invert) {
  sort(clades.begin(), clades.end(), [](const Clade& a, const Clade& b){ return a.size() < b.size(); });
  
  for (Clade& clade : clades){
    clade.score(scorer, clades, cladetaxa);
  }
  double score = clades.back().score(scorer, clades, cladetaxa);
  if (invert) { score = -score; }
  //BOOST_LOG_TRIVIAL(info) << "Score: " << format("%f") % score;
  cout << "Score: " << score;
  newick_tree = clades.back().newick_str(scorer, clades) ;
  cout << "Tree: " << newick_tree << endl;
   

}
