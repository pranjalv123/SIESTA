#include "TripartitionScorer.hpp"
#include "CladeExtractor.hpp"
#include <util/Logger.hpp>

#include <limits>
#include <fstream>
#include <cmath>
//#include <gperftools/profiler.h>

TripartitionScorerFactory::map_type* TripartitionScorerFactory::mymap;

double TripartitionScorer::get_score(clade_bitset& clade) {
  if(score_map.count(clade)){
    double outval = score_map[clade];
    return outval;
  }
  return nan("");
}

void TripartitionScorer::set_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2) {
  score_map[clade] = score;
  subclade_map.emplace(clade, make_pair(a1, a2));
}


void TripartitionScorer::add_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2) {
  DEBUG << "ADDING " << clade.str() << "\t" << a1.str() << "\t" << a2.str() << "\t" << score << endl;
  score_map[clade] = score;
  subclade_list_map[clade].emplace_back(make_pair(a1, a2));
}

void TripartitionScorer::clear_scores(clade_bitset& clade) {
  DEBUG << "CLEARING " << clade.str() << endl;
  subclade_list_map[clade].clear();
}


pair<clade_bitset, clade_bitset>& TripartitionScorer::get_subclades(clade_bitset& clade, vector<ScorableClade>& clades) {
  if(subclade_map.count(clade) == 0){
    Clade c(ts, clade);
    ERR << c.str() << " doesn't have subclades!" << endl;
    assert(false);
  }
  pair<clade_bitset, clade_bitset>& outval = subclade_map.at(clade);
  return outval;
}


vector<pair<clade_bitset, clade_bitset> >& TripartitionScorer::get_subclade_lists(clade_bitset& clade) {
  if(subclade_list_map.count(clade) == 0){
    Clade c(ts, clade);
    ERR << c.str() << " doesn't have subclades!" << endl;
    assert(false);
  }
  vector<pair<clade_bitset, clade_bitset> >& outval = subclade_list_map.at(clade);
  return outval;
}



double TripartitionScorer::adjust_final_score(double score) {
  return score; 
}

