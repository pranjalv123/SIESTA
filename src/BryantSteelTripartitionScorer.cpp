#include "TripartitionScorer.hpp"
#include "BryantSteelTripartitionScorer.hpp"
#include "CladeExtractor.hpp"
#include "Logger.hpp"

#include <limits>
#include <fstream>
#include <cmath>

DEF_SCORER(BryantSteelTripartitionScorer);

BryantSteelTripartitionScorer::BryantSteelTripartitionScorer(TaxonSet& ts) :
  TripartitionScorer(ts),
  qd(*QuartetDict::cl(ts))
{
  unordered_set<Clade>& clades = CladeExtractor::get_clades();
  for (const Clade& clade : clades) {
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
