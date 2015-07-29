#include "TripartitionScorer.hpp"

#include <boost/log/trivial.hpp>
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

double TripartitionScorer::get_score(clade_bitset& clade) {
  if(score_map.count(clade)){
    return score_map[clade];
  }
  return numeric_limits<double>::infinity();
}

void TripartitionScorer::set_score(clade_bitset& clade, double score, clade_bitset& a1, clade_bitset& a2) {
  score_map[clade] = score;
  //  BOOST_LOG_TRIVIAL(debug) << "SAVING" << clade << "\t" << (int)score << endl;
  subclade_map.emplace(clade, make_pair(a1, a2));

}

pair<clade_bitset, clade_bitset>& TripartitionScorer::get_subclades(clade_bitset& clade, vector<Clade>& clades) {
  if(subclade_map.count(clade) == 0){
    Clade c(ts, clade);
    BOOST_LOG_TRIVIAL(error) << c.str() << " doesn't have subclades!" << endl;
    //    c.score(*this, clades, cladetaxa);
    assert(false);
  }
  return subclade_map.at(clade);
}
