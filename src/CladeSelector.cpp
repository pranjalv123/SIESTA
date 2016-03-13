#include "CladeSelector.hpp"
#include "Logger.hpp"
#include <algorithm>
#include <iostream>
#include <vector>



double CladeSelector::run(bool invert) {  
 
  sort(clades.begin(), clades.end(), [](const Clade& a, const Clade& b){ return a.size() < b.size(); });

  INFO << "Scoring " << clades.size() << " clades" << endl;


  int current_size = 0;

  vector<vector<Clade*> > splitup(clades[clades.size()-1].size() + 1);

  for (Clade& c: clades) {
    splitup.at(c.size()).push_back(&c);
  }

  for (vector<Clade*> sublist : splitup) {
    INFO << "Processing " << sublist.size() << " clades of size " << current_size << endl;
#pragma omp parallel for
    for (size_t i = 0; i < sublist.size(); i++){
      Clade& clade = *(sublist[i]);
      DEBUG << clade.str() << endl;
      clade.score(scorer, clades, cladetaxa);
    }
    
    DEBUG << "Finished processing clades of size " << current_size << endl;
    current_size++;
  }

  
  
  double score = clades.back().score(scorer, clades, cladetaxa);
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
