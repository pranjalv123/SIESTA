#include "CladeSelector.hpp"
#include <algorithm>
#include "boost/boost/format.hpp"
#include "boost/boost/log/trivial.hpp"

using boost::format;
using boost::io::group;

void CladeSelector::run(bool invert) {
  sort(clades.begin(), clades.end(), [](Clade& a, Clade& b){ return a.size() < b.size(); });
  
  for (Clade& clade : clades){
    clade.score(scorer, clades, cladetaxa);
  }
  double score = clades.back().score(scorer, clades, cladetaxa);
  if (invert) { score = -score; }
  BOOST_LOG_TRIVIAL(info) << "Score: " << format("%f") % score;
  newick_tree = clades.back().newick_str(scorer, clades) ;
  BOOST_LOG_TRIVIAL(info) << "Tree: " << newick_tree << endl;
   

}
