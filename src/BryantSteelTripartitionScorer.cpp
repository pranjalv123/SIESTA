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

  if (Options::get("dominant")) {
    DEBUG << "using dominant quartet only" << endl;
    for (Taxon a = 0; a < ts.size(); a++) {
      for (Taxon b = 0; b < a; b++) {
	for (Taxon c = 0; c < b; c++) {
	  for (Taxon d = 0; d < c; d++) {
	    double q1 = qd(a,b,c,d);
	    double q2 = qd(a,c,b,d);
	    double q3 = qd(a,d,b,c);
	    double val = max(max(q1, q2), q3);
	    qd.set(a,b,c,d,0);
	    qd.set(a,c,b,d,0);
	    qd.set(a,d,b,c,0);
	    
	    if (q1 == val) {
	      qd.set(a,b,c,d,1);	      
	    }
	    
	    if (q2 == val) {
	      qd.set(a,c,b,d,1);	      
	    }
	    
	    if (q3 == val) {
	      qd.set(a,d,b,c,1);	      
	    }
	  }
	}
      }
    }
  }
  
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
	DEBUG << "Taxon pair" << i << "\t" << j << endl;
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

double BryantSteelTripartitionScorer::score(const Tripartition<ScorableClade>& t) {

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
