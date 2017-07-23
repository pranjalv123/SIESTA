#ifndef __SCORABLE_CLADE_HPP
#define __SCORABLE_CLADE_HPP

#include <vector>
#include "TripartitionScorer.hpp"
#include <Clade.hpp>
#include <util/Options.hpp>


typedef boost::multi_array<double, 2> twod_mat;

class TripartitionScorer;

class ScorableClade : public Clade {
  static int invert;
  static int initialized;
public:

  int myIndex;
  
  using Clade::Clade;
  
  ScorableClade(const Clade& c)
    : Clade (c),
      myIndex(-1)
  {}
  
  string newick_str(vector<ScorableClade>& clades, vector<double>& supports, size_t index);
  string newick_str(TripartitionScorer& scorer, vector<ScorableClade>& clades);
  vector<string> all_newick_strs(TripartitionScorer& scorer, vector<ScorableClade>& clades, twod_mat& mat, unordered_map<clade_bitset, vector<string> >& cache);

  __int128 optimal_subtree_count(TripartitionScorer& scorer, unordered_map<clade_bitset, __int128 >& cache);
  double defective_subtree_count(TripartitionScorer& scorer, twod_mat& mat, vector<ScorableClade>& clades, unordered_map<clade_bitset, int>& clade_indices, double defect, unordered_map< clade_bitset, unordered_map<double, double > >& cache);
  __int128 appearances_in_optimal_trees(TripartitionScorer& scorer, unordered_map<clade_bitset, __int128 >& count_cache) const;
  double appearances_in_defective_trees(TripartitionScorer& scorer, double defect, unordered_map<clade_bitset, unordered_map<double, double> >& defective_cache);

  double score(TripartitionScorer& scorer, vector<ScorableClade>& clades, unordered_map<clade_bitset, int>& cladetaxa, twod_mat* mat);
};


struct ScorableTripartition {
  ScorableClade a1, a2, rest;
  ScorableTripartition(const TaxonSet& ts, const ScorableClade& clade, const ScorableClade& subclade);
  string str() const;
};

#endif
