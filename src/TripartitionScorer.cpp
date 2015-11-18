#include "TripartitionScorer.hpp"
#include "CladeExtractor.hpp"

#include <limits>
#include <fstream>
#include <gperftools/profiler.h>

BryantSteelTripartitionScorer::BryantSteelTripartitionScorer(TaxonSet& ts, QuartetDict& qd, vector<Clade>& clades) :
  TripartitionScorer(ts),
  qd(qd)
{
  for (Clade& clade : clades) {
    vector<Taxon> nonmembers;
    for(size_t i = 0; i < ts.size(); i++) {
      if (!clade.contains(i))
  	nonmembers.push_back(i);
    }

    map<pair<Taxon, Taxon>, double>& mp = W[clade.get_taxa()];

    if (clade.size() < 2) {
      continue;
    }
    
    for (size_t i = 0; i < nonmembers.size(); i++) {
      for (size_t j = i+1; j < nonmembers.size(); j++) {

  	double d = 0;
  	for (Taxon k : clade) {
  	  for (Taxon l : clade) {
  	    if (k > l)
  	      d += qd(nonmembers[i],nonmembers[j],k,l);
  	  }
  	}
	
  	mp[make_pair(nonmembers[i], nonmembers[j])] = d;
  	mp[make_pair(nonmembers[j], nonmembers[i])] = d;
      }
    }
  }
}

double BryantSteelTripartitionScorer::score(const Tripartition& t) {

  double val = 0;


  for (Taxon c : t.a2) 
    for (Taxon d : t.rest) 
      val +=  W[t.a1.get_taxa()][make_pair(c,d)];
    


  for (Taxon c : t.a1) 
    for (Taxon d : t.rest) 
      val +=  W[t.a2.get_taxa()][make_pair(c,d)];
  

  for (Taxon c : t.a2) 
    for (Taxon d : t.a2)
      if (c > d)
	val += W[t.a1.get_taxa()][make_pair(c, d)];

  return val;
}

double DPTripartitionScorer::score(const Tripartition& t) {
  double val = 0;
  
  for(Taxon a : t.a1) 
    for (Taxon b: t.a2) 
      for (Taxon c : t.rest) 
	for (Taxon d : t.rest) 
	  if (c > d)
	    val += qd(a,b,c,d);


  for(Taxon a : t.a1) 
    for (Taxon b: t.a1) 
      if (a > b) 
	for (Taxon c : t.a2) 
	  for (Taxon d : t.a2) 
	    if (c > d)
	      val -= qd(a,b,c,d);

  return val;
}


RFTripartitionScorer::RFTripartitionScorer(TaxonSet& ts, string treesfile) :
  TripartitionScorer(ts)
{
  string tree;
  ifstream file(treesfile);
  while(getline(file, tree)) {
    unordered_set<Taxon> tree_taxa;
    unordered_set<Clade> clades = CladeExtractor::extract(ts, tree, tree_taxa);
    unordered_set<Clade> clade_complements;
    Clade tree_clade(ts, tree_taxa);
    for (const Clade& clade : clades) {
      Clade comp(tree_clade.minus(clade));
      clade_weights[Bipartition(clade, comp)] += 1;
    }    
  }


}

RFTripartitionScorer::RFTripartitionScorer(TaxonSet& ts, vector<string> trees) :
  TripartitionScorer(ts)
{
  for (string& tree : trees) {
    unordered_set<Taxon> tree_taxa;
    unordered_set<Clade> clades = CladeExtractor::extract(ts, tree, tree_taxa);
    unordered_set<Clade> clade_complements;
    Clade tree_clade(ts, tree_taxa);
    for (const Clade& clade : clades) {
      Clade comp(tree_clade.minus(clade));
      clade_weights[Bipartition(clade, comp)] += 1;
    }
  }
}

bool RFTripartitionScorer::matches(const Tripartition& t, const Bipartition& bp) {
  if (t.a1.contains(bp.a1) && t.a1.complement().contains(bp.a2) && !t.rest.contains(bp.a2))
    return true;

  if (t.a1.contains(bp.a2) && t.a1.complement().contains(bp.a1) && !t.rest.contains(bp.a1))
    return true;

  if (t.a2.contains(bp.a1) && t.a2.complement().contains(bp.a2) && !t.rest.contains(bp.a2))
    return true;

  if (t.a2.contains(bp.a2) && t.a2.complement().contains(bp.a1) && !t.rest.contains(bp.a1))
    return true;

  return false;
}

double RFTripartitionScorer::score(const Tripartition& t) {
  double weight = 0;
  for (auto i: clade_weights) {
    const Bipartition& bp = i.first;
    double c_weight = i.second;
    
    if (matches(t, bp)) {
      weight += c_weight;
    }
  }
  return weight;
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
    //    BOOST_LOG_TRIVIAL(error) << c.str() << " doesn't have subclades!" << endl;
    //    c.score(*this, clades, cladetaxa);
    assert(false);
  }
  return subclade_map.at(clade);
}

