#include "TripartitionScorer.hpp"
#include "FastRFTripartitionScorer.hpp"
#include "CladeExtractor.hpp"
#include "Logger.hpp"

#include <limits>
#include <fstream>
#include <cmath>
DEF_SCORER(FastRFTripartitionScorer);



FastRFTripartitionScorer::FastRFTripartitionScorer(TaxonSet& ts) :
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
  INFO << "Total Weight: " << total_weight << endl;

  for (auto& i : clade_weights) {
    bipartitions.push_back(&(i.first));
    weights.push_back(i.second);
  }
  
  unordered_set<Clade>& clades = CladeExtractor::get_clades();

  INFO << "Making possibilities lists\n";


  int n = 0;
  for (const Clade& clade : clades) {
    BitVectorFixed& bvf_a1 = (*(possibleA1.emplace(clade, weights.size()).first)).second;
    BitVectorFixed& bvf_a2 = (*(possibleA2.emplace(clade, weights.size()).first)).second;
    BitVectorFixed& bvf_b1 = (*(possibleB1.emplace(clade, weights.size()).first)).second;
    BitVectorFixed& bvf_b2 = (*(possibleB2.emplace(clade, weights.size()).first)).second;
    
    for (unsigned int i = 0; i < bipartitions.size(); i++) {
      const Bipartition& bp = *(bipartitions[i]);
      if (clade.contains(bp.a1) && clade.overlap(bp.a2).size() == 0) {
	bvf_a1.set(i);
      }
      if (clade.overlap(bp.a1).size() > 0) {
	bvf_b1.set(i);
      }
      if (clade.contains(bp.a2) && clade.overlap(bp.a1).size() == 0) {
	bvf_a2.set(i);
      }
      if (clade.overlap(bp.a2).size() > 0) {
	bvf_b2.set(i);
      }
    }
    n++;
    if (n%10000 == 0) {
      INFO << "prepared " << n << "/" << clades.size() <<  endl;
    }
    
  }
}

int FastRFTripartitionScorer::addSourceTree(string tree) {
  unordered_set<Taxon> tree_taxa;
  unordered_set<Clade> clades = CladeExtractor::extract(ts, tree, tree_taxa);
  unordered_set<Clade> clade_complements;
  Clade tree_clade(ts, tree_taxa);
  int n = 0;
  
  for (const Clade& clade : clades) {
    Clade comp(tree_clade.minus(clade));
    //    Clade comp(clade.complement());
    if (comp.size() && clade.size()) {
      clade_weights[Bipartition(clade, comp)] += 1;
      n++;
    }
  }

  return n;
}

bool FastRFTripartitionScorer::matches(const Tripartition<ScorableClade>& t, const Bipartition& bp) {

  if (t.a1.contains(bp.a1) && (t.a1.overlap(bp.a2).size() == 0) && (t.a2.overlap(bp.a2).size() > 0)) {
    //    DEBUG << "match 1\n";
    return true;
  }
  if (t.a2.contains(bp.a1) && (t.a2.overlap(bp.a2).size() == 0) && (t.a1.overlap(bp.a2).size() > 0)) {
    //    DEBUG << "match 2\n";
    return true;
  }

  if (t.a1.contains(bp.a2) && (t.a1.overlap(bp.a1).size() == 0) && (t.a2.overlap(bp.a1).size() > 0)) {
    //    DEBUG << "match 3\n";
    return true;
  }
  if (t.a2.contains(bp.a2) && (t.a2.overlap(bp.a1).size() == 0) && (t.a1.overlap(bp.a1).size() > 0)) {
    //    DEBUG << "match 4\n";
    return true;
  }
  
  
  return false;
}

double FastRFTripartitionScorer::score(const Tripartition<ScorableClade>& t) {
  double weight = 0;

  BitVectorFixed& t1_a1 = possibleA1.at(t.a1);
  BitVectorFixed& t1_a2 = possibleA2.at(t.a1);
  BitVectorFixed& t1_b1 = possibleB1.at(t.a1);
  BitVectorFixed& t1_b2 = possibleB2.at(t.a1);

  BitVectorFixed& t2_a1 = possibleA1.at(t.a2);
  BitVectorFixed& t2_a2 = possibleA2.at(t.a2);
  BitVectorFixed& t2_b1 = possibleB1.at(t.a2);
  BitVectorFixed& t2_b2 = possibleB2.at(t.a2);

  BitVectorFixed combinations =
    (t1_a1 & t2_b2) | (t2_a1 & t1_b2) | (t1_a2 & t2_b1) | (t2_a2 & t1_b1);
  
  for (int i : combinations) {
    weight += weights[i];
  }
  
  
  return weight;
}

double FastRFTripartitionScorer::adjust_final_score(double score) {
  PROGRESS << "Raw Score: " << score << endl;
  return ((total_weight - score) ) * 2;
}
