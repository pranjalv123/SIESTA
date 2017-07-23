#include "CladeSelector.hpp"
#include "ScorableClade.hpp"
#include <util/Logger.hpp>
#include "Options.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
#ifdef ENABLE_PROFILING
#include <gperftools/profiler.h>
#endif



double BasicCladeSelector::run(bool invert, twod_mat* mat) {
#ifdef ENABLE_PROFILING
  bool profile = Options::get("profile");
#endif

  unordered_map<clade_bitset, int> clade_indices;
  
  INFO << "Sorting " << clades.size() << " clades" << endl;

  sort(clades.begin(), clades.end(), [](const ScorableClade& a, const ScorableClade& b){ return a.size() < b.size(); });
  
  
  for ( size_t i = 0; i < clades.size(); i++ ) {
    clades[i].myIndex = i;
    clade_indices[clades[i].taxa] = i;
  }
  
  INFO << "Scoring " << clades.size() << " clades" << endl;


  int current_size = 0;

  vector<vector<ScorableClade*> > splitup(clades[clades.size()-1].size() + 1);
  
  for (ScorableClade& c: clades) {
    splitup.at(c.size()).push_back(&c);
  }

  for (vector<ScorableClade*> sublist : splitup) {
    INFO << "Processing " << sublist.size() << " clades of size " << current_size << endl;

    for (size_t i = 0; i < sublist.size(); i++) {
      ScorableClade& clade = *(sublist[i]);
      //      DEBUG << clade.str() << endl;
      clade.score(scorer, clades, clade_indices, mat);
    }
    #ifdef ENABLE_PROFILING
    if (profile)
      ProfilerFlush();
    #endif

    DEBUG << "Finished processing clades of size " << current_size << endl;
    current_size++;
  }

  
  
  double score = clades.back().score(scorer, clades, clade_indices, mat);
  if (invert) { score = -score; }
  //BOOST_LOG_TRIVIAL(info) << "Score: " << format("%f") % score;
  cout << "Score: " << scorer.adjust_final_score(score) << endl;
  newick_tree = clades.back().newick_str(scorer, clades) ;
  cout << "Tree: " << newick_tree << endl;
  
  return score;
}







// double MPICladeSelector::run(bool invert) {  
  
//   sort(clades.begin(), clades.end(), [](const Clade& a, const Clade& b){ return a.size() < b.size(); });

//   INFO << "Scoring " << clades.size() << " clades" << endl;

//   int n = 0;
  
//   for (Clade& clade : clades){
//     DEBUG << clade.str() << endl;
//     clade.score(scorer, clades, cladetaxa);
//     if (n % 1000 == 0) {
//       INFO << "Scored " << n << "/" << clades.size() << endl;
//     }
//     n++;
//   }
//   double score = clades.back().score(scorer, clades, cladetaxa);
//   if (invert) { score = -score; }
//   //BOOST_LOG_TRIVIAL(info) << "Score: " << format("%f") % score;
//   cout << "Score: " << scorer.adjust_final_score(score) << endl;
//   newick_tree = clades.back().newick_str(scorer, clades) ;
//   cout << "Tree: " << newick_tree << endl;
  
//   return score;
// }
