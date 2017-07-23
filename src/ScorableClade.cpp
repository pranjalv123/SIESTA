#include "ScorableClade.hpp"

#include<algorithm>
#include<limits>

int ScorableClade::invert = 1;
int ScorableClade::initialized = 0;

string ScorableClade::newick_str(vector<ScorableClade>& clades, vector<double>& supports, size_t index) {
  if (size() == 0) {
    return "";
  }
  if (size() == 1) {
    return ts[*begin()];
  }
  if (size() == 2) {
    stringstream ss;

    vector<Taxon> tv;
    for (Taxon t : *this) {
      tv.push_back(t);
    }
    
    ss << "("<<ts[tv[0]] << "," << ts[tv[1]] << "):" << supports[index] ;    
    
    return ss.str();
  }

  ScorableClade remaining(*this);
  
  vector<size_t> indices;

  auto it = clades.begin() + index;
  DEBUG << str() << endl;
  while (remaining.size()) {
    it++;
    
    it = find_if(it, clades.end(), [ & ] (ScorableClade& c){return remaining.contains(c) && !(c == *this);});
    
    ScorableClade& sc = *it;
    DEBUG << endl;
    DEBUG << remaining.str() << remaining.taxa.str() << endl;
    DEBUG << sc.str() << " " << sc.taxa.str() << endl;
    DEBUG << remaining.minus(sc).str() << " " << remaining.minus(sc).taxa.str() << endl;

    remaining = remaining.minus(sc);
    DEBUG << remaining.str() << endl;
    DEBUG << distance( clades.begin(), it) << endl;
    indices.push_back(distance( clades.begin(), it));
    
  }

  
  stringstream ss;
  ss << "(";

  bool started = false;

  
  DEBUG << "Clade " << str() << " has subclades:" << endl;
  for (size_t ix : indices) {
    DEBUG << ix << endl;
    DEBUG << clades[ix].str() << endl;
  }

  if (indices.size() > 2) {
    INFO << "Node with degree " << indices.size() << endl;
  }
  
  for (size_t ix : indices) {
    
    if (started)
      ss << ",";
    ss << clades[ix].newick_str(clades, supports, ix);
    started = true;
  }
  ss << "):" << supports[index];
    
  //  BOOST_LOG_TRIVIAL(debug) << str() << c1.str() << c2.str() << (int)scorer.get_score(taxa) <<endl;

  return ss.str();
}

string ScorableClade::newick_str(TripartitionScorer& scorer, vector<ScorableClade>& clades) {
  //  BOOST_LOG_TRIVIAL(info) << str() << endl;
  if (size() == 0) {
    return "";
  }
  if (size() == 1) {
    return ts[*begin()];
  }
  if (size() == 2) {
    stringstream ss;

    vector<Taxon> tv;
    for (Taxon t : *this) {
      tv.push_back(t);
    }
    
    ss << "("<<ts[tv[0]] << "," << ts[tv[1]] << ")" ;    
    
    return ss.str();

  }


  stringstream ss;
  pair<clade_bitset, clade_bitset>& subclades = scorer.get_subclades(taxa, clades);

  ScorableClade c1(ts, subclades.first);
  ScorableClade c2(ts, subclades.second);
  
  //  BOOST_LOG_TRIVIAL(debug) << str() << c1.str() << c2.str() << (int)scorer.get_score(taxa) <<endl;

  Tripartition<ScorableClade> tp(ts, *this, c1);
  
  ss << "(" << c1.newick_str(scorer, clades) << "," << c2.newick_str(scorer, clades) << ")";
  return ss.str();
}


double ScorableClade::defective_subtree_count(TripartitionScorer& scorer, twod_mat& mat, vector<ScorableClade>& clades, unordered_map<clade_bitset, int>& clade_indices, double defect, unordered_map< clade_bitset, unordered_map<double, double> >& cache) {

  int index = clade_indices[taxa];

  if (size() <= 2) {
    if (defect == 0) {
      cache[taxa][defect] = 1;
      return 1;
    }
    cache[taxa][defect] = 0;
    return 0;
  }
  
  if (cache.count(taxa) && cache[taxa].count(defect)) {
    return cache[taxa][defect];
  }

  double score = scorer.get_score(taxa);

  double count = 0;

  
  for (int i = 0; i < clades.size(); i++) {
    double val = mat[index][clade_indices[clades[i].taxa]];
    if (isnan(val))
      continue;
    
    Tripartition<ScorableClade> tp(ts, *this, clades[i]);
        

    if (tp.a1.taxa.ffs() < tp.a2.taxa.ffs())
      continue;
    
    
    if (abs(val - score) <= defect) {
      double residual = defect - abs(val - score);
      
      for (double v = 0; v <= residual; v++)
	count +=
	  tp.a1.defective_subtree_count(scorer, mat, clades, clade_indices, v, cache) * tp.a2.defective_subtree_count(scorer, mat, clades, clade_indices, residual - v, cache) ;
    }    
  }
  
  cache[taxa][defect] = count;
  return count;
}


__int128 ScorableClade::optimal_subtree_count(TripartitionScorer& scorer, unordered_map<clade_bitset, __int128 >& cache) {
  
  if (cache.count(taxa)) {
    return cache[taxa];
  }

  if (size() <= 2) {
    cache[taxa] = 1;
    return 1;
  }

  
  vector<pair<clade_bitset, clade_bitset> >& subclades_list = scorer.get_subclade_lists(taxa);

  __int128 count = 0;
  
  for (auto& subclades : subclades_list) {

    ScorableClade c1(ts, subclades.first);
    ScorableClade c2(ts, subclades.second);
    DEBUG << size() << "\t" << c1.size() << "\t" << c2.size() << endl;
    assert(c1.size() < size());
    assert(c2.size() < size());
    assert(c1.size() + c2.size() == size());
    count += c1.optimal_subtree_count(scorer, cache) * c2.optimal_subtree_count(scorer, cache);
    
  }
  cache[taxa] = count;
  return count;
}

__int128 ScorableClade::appearances_in_optimal_trees(TripartitionScorer& scorer, unordered_map<clade_bitset, __int128 >& count_cache) const { 

  return count_cache[taxa] * count_cache[complement().taxa];
}

double ScorableClade::appearances_in_defective_trees(TripartitionScorer& scorer, double defect, unordered_map<clade_bitset, unordered_map<double, double> >& defective_cache) {
  
  double count = 0;
  for (int i = 0; i <= defect; i++) {
    for (int j = 0; j <= defect - i; j++) {
      count += defective_cache[taxa][i] * defective_cache[complement().taxa][j];
    }
  }
  return count;
}

vector<string> ScorableClade::all_newick_strs(TripartitionScorer& scorer, vector<ScorableClade>& clades, twod_mat& mat, unordered_map<clade_bitset, vector<string> >& cache) {

  if (cache.count(taxa)) {
    return cache[taxa];
  }
  
  vector<string>& output = cache[taxa];
  //  BOOST_LOG_TRIVIAL(info) << str() << endl;
  if (size() == 0) {
    output.push_back("");
    cache[taxa] = output;
    return output;
  }
  if (size() == 1) {
    output.push_back(ts[*begin()]);
    cache[taxa] = output;
    return output;
  }
  if (size() == 2) {

    vector<Taxon> tv;
    for (Taxon t : *this) {
      tv.push_back(t);
    }
    stringstream ss;    
    ss << "("<<ts[tv[0]] << "," << ts[tv[1]] << ")" ;    

    output.push_back(ss.str());
    cache[taxa] = output;
    return output;

  }




  vector<pair<clade_bitset, clade_bitset> >& subclades_list = scorer.get_subclade_lists(taxa);

  if (size() == ts.size()) {
    subclades_list.erase(subclades_list.begin() + 1, subclades_list.end());
  }
  
  for (auto& subclades : subclades_list) {

    ScorableClade c1(ts, subclades.first);
    ScorableClade c2(ts, subclades.second);

    DEBUG << "Subclades of " << str() << " "  << c1.str() << " "<< c2.str() << endl;
    
    //  BOOST_LOG_TRIVIAL(debug) << str() << c1.str() << c2.str() << (int)scorer.get_score(taxa) <<endl;

    for (string s1 : c1.all_newick_strs(scorer, clades, mat, cache)) {
      for (string s2 : c2.all_newick_strs(scorer, clades, mat, cache)) {
	stringstream ss;
	ss << "(" << s1 << "," << s2 << ")";
	output.push_back(ss.str());
	DEBUG << ss.str() << endl;
      }
    }
  }
  cache[taxa] = output;
  return output;
}




double ScorableClade::score(TripartitionScorer& scorer, vector<ScorableClade>& clades, unordered_map<clade_bitset, int>& clade_indices, twod_mat* mat) {
  double value;

  
  if (size() == 1) {
    value = 0;
    Clade eclade (ts);
    scorer.set_score(taxa, value, taxa, eclade.taxa);
    return value;
  }

  value = scorer.get_score(taxa);
  if (!std::isnan(value)) {
    return value;
  }
  clade_bitset sub1(ts.size()), sub2(ts.size());

  if (!initialized) {
    if (Options::get("maximize")) {
      invert = -1;
    }
    else {
      invert = 1;
    }
    initialized = 1;
  }
  
  if (size() == 2) {
    ScorableClade c1(ts);
    c1.add(*begin());
    Tripartition<ScorableClade> tp(ts, *this, c1);
    value = invert * scorer.score(tp);

    scorer.set_score(taxa, value, tp.a1.taxa, tp.a2.taxa);
  }
  else {
    vector<double> values(clades.size());
    DEBUG << "scoring: " << taxa.str() << endl;
    for (size_t i = 0; i < clades.size(); i++) {
      values[i] = (double)nan("");
      ScorableClade& subclade = clades[i];
      if (subclade.size() >= size() || !contains(subclade) || subclade.size() == 0 )
	continue;

      Tripartition<ScorableClade> tp(ts, *this, subclade);

      if (clade_indices.count(tp.a1.taxa) == 0 || clade_indices.count(tp.a2.taxa) == 0)
	continue;

      double score = invert * scorer.score(tp) + tp.a1.score(scorer, clades, clade_indices, mat) + tp.a2.score(scorer, clades, clade_indices, mat);
      values[i] = score;
      if (std::isnan(value) || (score < value) ) {
	
	
	value = score;
	sub1 = tp.a1.taxa;
	sub2 = tp.a2.taxa;
	
	
	DEBUG << "value: " << fixed << value << "\t" << score << "\t" << values[i] << endl;
      }
      if (mat) {
	(*mat)[clade_indices[taxa]][clade_indices[tp.a1.taxa]] = score;
	(*mat)[clade_indices[taxa]][clade_indices[tp.a2.taxa]] = score;
      }
    }
    for (int i = 0; i < clade_indices[taxa]; i++) {
      if (values[i] == value) {
	Tripartition<ScorableClade> tp(ts, *this, clades[i]);

	if (tp.a1.taxa.ffs() > tp.a2.taxa.ffs()) {
	  DEBUG << str() << clades[i].str() << endl;
	  DEBUG << str() << " "  << tp.a1.str() << " " << tp.a2.str() << endl;	  
	  assert(tp.a1.size() + tp.a2.size() == size());
	  scorer.add_score(taxa, value, tp.a1.taxa, tp.a2.taxa);
	}
      }       
    }
    DEBUG << scorer.get_subclade_lists(taxa).size() << endl;
    scorer.set_score(taxa, value, sub1, sub2);    
  }
  
  
  return value;
  
}


ScorableTripartition::ScorableTripartition(const TaxonSet& ts, const ScorableClade& clade, const ScorableClade& subclade) :
  a1(clade.minus(subclade)),
  a2(subclade),
  rest(clade.complement())
{}


string ScorableTripartition::str()  const {  
  assert(a1.overlap(rest).size() == 0);
  assert(a2.overlap(rest).size() == 0);
  assert(a2.overlap(a1).size() == 0);
  return "{" + a1.str() + "/" + a2.str() + "/" + rest.str() + "}";
}
