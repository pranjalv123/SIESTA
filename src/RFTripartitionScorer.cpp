#include "TripartitionScorer.hpp"
#include "RFTripartitionScorer.hpp"
#include "CladeExtractor.hpp"
#include "Logger.hpp"

#include <limits>
#include <fstream>
#include <cmath>
DEF_SCORER(RFTripartitionScorer);



RFTripartitionScorer::RFTripartitionScorer(TaxonSet& ts) :
  TripartitionScorer(ts)
{
  string treesfile;
  Options::get("g genetrees", &treesfile);
  string tree;
  ifstream file(treesfile);
  total_weight = 0;
  n_trees = 0;
  while(getline(file, tree)) {
    total_weight += addSourceTree(tree);
    n_trees ++;
    DEBUG << n_trees << "\t" << total_weight << endl;
  }
  for (auto& i: clade_weights) {
    DEBUG << i.first.str() << "\t" << i.second << endl;
  }
  cout << "Total Weight: " << total_weight << endl;
}

int RFTripartitionScorer::addSourceTree(string tree) {
  unordered_set<Taxon> tree_taxa;
  unordered_set<Clade> clades = CladeExtractor::extract(ts, tree, tree_taxa);
  unordered_set<Clade> clade_complements;
  Clade tree_clade(ts, tree_taxa);
  for (const Clade& clade : clades) {
    Clade comp(tree_clade.minus(clade));
    //    Clade comp(clade.complement());
    if (comp.size() && clade.size())
      clade_weights[Bipartition(clade, comp)] += 1;
  }

  return clades.size();
}

bool RFTripartitionScorer::matches(const Tripartition& t, const Bipartition& bp) {
  // if (t.a1.contains(bp.a1) && t.a1.complement().contains(bp.a2) && !t.rest.contains(bp.a2))
  //   return true;

  // if (t.a1.contains(bp.a2) && t.a1.complement().contains(bp.a1) && !t.rest.contains(bp.a1))
  //   return true;

  // if (t.a2.contains(bp.a1) && t.a2.complement().contains(bp.a2) && !t.rest.contains(bp.a2))
  //   return true;

  // if (t.a2.contains(bp.a2) && t.a2.complement().contains(bp.a1) && !t.rest.contains(bp.a1))
  //   return true;

  if (t.a1.contains(bp.a1) && (t.a1.overlap(bp.a2).size() == 0) && (t.a2.overlap(bp.a2).size() > 0)) {
    DEBUG << "match 1\n";
    return true;
  }
  if (t.a2.contains(bp.a1) && (t.a2.overlap(bp.a2).size() == 0) && (t.a1.overlap(bp.a2).size() > 0)) {
    DEBUG << "match 2\n";
    return true;
  }

  if (t.a1.contains(bp.a2) && (t.a1.overlap(bp.a1).size() == 0) && (t.a2.overlap(bp.a1).size() > 0)) {
    DEBUG << "match 3\n";
    return true;
  }
  if (t.a2.contains(bp.a2) && (t.a2.overlap(bp.a1).size() == 0) && (t.a1.overlap(bp.a1).size() > 0)) {
    DEBUG << "match 4\n";
    return true;
  }
  
  
  return false;
}

double RFTripartitionScorer::score(const Tripartition& t) {
  double weight = 0;
  for (auto i: clade_weights) {
    const Bipartition& bp = i.first;
    double c_weight = i.second;
    
    if (matches(t, bp)) {
      DEBUG << t.str() << " matches " << bp.str() << endl;
      weight += c_weight;
    }
  }
  DEBUG << "weight " << t.str() << " " << weight << endl;
  return weight;
}

double RFTripartitionScorer::adjust_final_score(double score) {
  PROGRESS << "Raw Score: " << score << endl;
  return ((total_weight - score) - n_trees) * 2;
}
