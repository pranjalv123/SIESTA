#include "TripartitionScorer.hpp"
#include "AstralTripartitionScorer.hpp"
#include "CladeExtractor.hpp"
#include "Logger.hpp"

#include <newick.hpp>
#include <BitVector.hpp>
#include <limits>
#include <fstream>
#include <cmath>

DEF_SCORER(ASTRALTripartitionScorer);

ASTRALTripartitionScorer::ASTRALTripartitionScorer(TaxonSet& ts) :
  TripartitionScorer(ts)
{
  string treesfile;
  Options::get("g genetrees", &treesfile);

  ifstream file(treesfile);
  string tree;
  
  while(getline(file, tree)) {
    vector<Taxon> thistree;
    newick_to_postorder(tree, ts, thistree);       
    postorders.push_back(thistree);
    gttaxa.push_back(newick_to_taxa(tree,ts));

  }

}


double F(double a, double b, double c) {
  return a*b*c*(a+b+c-3)/2;
}


double qi(array<double, 3> c1, array<double, 3> c2, array<double, 3> c3) {
  return
    F(c1[0], c2[1], c3[2]) +
    F(c1[0], c2[2], c3[1]) +
    F(c1[1], c2[0], c3[2]) +
    F(c1[1], c2[2], c3[0]) +
    F(c1[2], c2[1], c3[0]) +
    F(c1[2], c2[0], c3[1]);    
}

double ASTRALTripartitionScorer::score(const Tripartition<ScorableClade>& t) {

  double val = 0;
  for (int ix = 0; ix < postorders.size(); ix++) {
    auto& po =  postorders[ix];
    vector<array<double, 3> > stack;
    for (Taxon u : po) {
      if (u >= 0) {
	stack.push_back({t.a1.contains(u) > 0, t.a2.contains(u) > 0, t.rest.contains(u)> 0} );	
      } else {
	if (u == -2) {
	  assert(stack.size());
	  array<double, 3> c1 = stack.back();
	  stack.pop_back();
	  assert(stack.size());
	  array<double, 3> c2 = stack.back();
	  stack.pop_back();

	  double x=  c1[0] + c2[0];
	  double y=  c1[1] + c2[1];
	  double z=  c1[2] + c2[2];
	  
	  if (gttaxa[ix].size() < ts.size())  {
	    
	    array<double, 3> c3 = {t.a1.taxa.overlap_size(gttaxa[ix].taxa) - x, t.a2.taxa.overlap_size(gttaxa[ix].taxa) - y, t.rest.taxa.overlap_size(gttaxa[ix].taxa) - z};
	    val += qi(c1, c2, c3);
	  }
	  else {
	    array<double, 3> c3 = {t.a1.size() - x, t.a2.size() - y, t.rest.size() - z};

	    val += qi(c1, c2, c3);
	  }	 
	  stack.push_back({x,y,z} );
	}
	else {
	  vector<array<double, 3> > children;

	  for (int i = 0; i < -u; i++) {	  
	    assert(stack.size());
	    children.push_back(stack.back());
	    stack.pop_back();
	  }
	  
	  double x = 0;
	  double y = 0;
	  double z = 0;
	
	  for (int i = 0; i < -u; i++) {
	    x += children[i][0];
	    y += children[i][1];
	    z += children[i][2];
	  }
	
	  array<double, 3> c3 = {t.a1.taxa.overlap_size(gttaxa[ix].taxa) - x, t.a2.taxa.overlap_size(gttaxa[ix].taxa) - y, t.rest.taxa.overlap_size(gttaxa[ix].taxa) - z};

	  children.push_back(c3);

	  for (int i = 0; i < children.size(); i++) {
	    for (int j = i+1; j < children.size(); j++) {
	      for (int k = j+1; k < children.size(); k++) {	      
		val += qi(children[i], children[j], children[k]);  
	      }
	    }
	  }
	  
	  stack.push_back({x,y,z} );	
	}
	
      }
    }
  }
  
  return val;
}

double ASTRALTripartitionScorer::adjust_final_score(double score) {
  return score/2;
}
