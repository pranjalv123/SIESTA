#include "TripartitionScorer.hpp"
#include "CladeExtractor.hpp"
#include "Logger.hpp"

#include <limits>
#include <fstream>
#include <cmath>
//#include <gperftools/profiler.h>

TripartitionScorerFactory::map_type* TripartitionScorerFactory::mymap;

double TripartitionScorer::get_score(clade_bitset& clade) {
  if(score_map.count(clade)){
    return score_map[clade];
  }
  return nan("");
}

void TripartitionScorer::set_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2) {
  score_map[clade] = score;
  subclade_map.emplace(clade, make_pair(a1, a2));

}

pair<clade_bitset, clade_bitset>& TripartitionScorer::get_subclades(clade_bitset& clade, vector<Clade>& clades) {
  if(subclade_map.count(clade) == 0){
    Clade c(ts, clade);
    ERR << c.str() << " doesn't have subclades!" << endl;
    assert(false);
  }
  return subclade_map.at(clade);
}


double TripartitionScorer::adjust_final_score(double score) {
  return score; 
}

