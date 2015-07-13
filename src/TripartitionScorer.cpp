#include "TripartitionScorer.hpp"
#include <limits>

double DPTripartitionScorer::score(const Tripartition& t) {
  double val = 0;
  
  for(Taxon a : t.a1.taxa_list) 
    for (Taxon b: t.a2.taxa_list) 
      for (Taxon c : t.rest.taxa_list) 
	for (Taxon d : t.rest.taxa_list) 
	  if (c > d)
	    val += qd(a,b,c,d);


  for(Taxon a : t.a1.taxa_list) 
    for (Taxon b: t.a1.taxa_list) 
      if (a > b) 
	for (Taxon c : t.a2.taxa_list) 
	  for (Taxon d : t.a2.taxa_list) 
	    if (c > d)
	      val -= qd(a,b,c,d);
  return val;
}
int misscount = 0;
double TripartitionScorer::get_score(bitset<128>& clade) {
  if(score_map.count(clade)){
    return score_map[clade];
  }
  // if (misscount % 1000 == 0) {
  //   cout << misscount <<  "\t" << score_map.size() << endl;
  // }
  misscount++;
  //  cout << Clade(ts, clade).str() << endl;
  return numeric_limits<double>::infinity();
}

void TripartitionScorer::set_score(bitset<128>& clade, double score, Clade& a1, Clade& a2) {
  
  score_map[clade] = score;
  subclade_map.emplace(clade, make_pair(a1, a2));
}

pair<Clade, Clade>& TripartitionScorer::get_subclades(bitset<128>& clade, vector<Clade>& clades) {
  if(subclade_map.count(clade) == 0){
    Clade c(ts, clade);
    cout << c.str() << " doesn't have subclades!" << endl;
    //    c.score(*this, clades, cladetaxa);
    assert(false);
  }
  return subclade_map.at(clade);
}
