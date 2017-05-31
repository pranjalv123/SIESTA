#include "TripartitionScorer.hpp"
#include "DPTripartitionScorer.hpp"
#include "CladeExtractor.hpp"
#include "Logger.hpp"

#include <limits>
#include <fstream>
#include <cmath>


DEF_SCORER(DPTripartitionScorer);


//make a quartet dictionary based on the command line arguments

DPTripartitionScorer::DPTripartitionScorer(TaxonSet& ts) :
  TripartitionScorer(ts),
  qd(*QuartetDict::cl(ts)) {}


//find quartet support of a tripartition
double DPTripartitionScorer::score(const Tripartition<ScorableClade>& t) {
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
