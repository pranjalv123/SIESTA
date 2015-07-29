#include "Clade.hpp"
#include "TripartitionScorer.hpp"
#include <sstream>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <boost/log/trivial.hpp>

Clade::Clade(TaxonSet& ts_, string& str) :
  ts(ts_)
{
  int i = 1;

  char* cladestr = &(str[1]);
  char* token;
  char* saveptr;

  while(token = strtok_r(cladestr, ",} ", &saveptr)) {
    cladestr = NULL;
    add(ts[string(&(token[0]))]);
  }
}

Clade::Clade(TaxonSet& ts_) :
  ts(ts_)
{}

Clade::Clade(const Clade& other) :
  taxa_list(other.taxa_list),
  taxa(other.taxa),
  ts(other.ts) {}

Clade::Clade(TaxonSet& ts_, clade_bitset taxa) :
  ts(ts_),
  taxa(taxa)
{
  for (int i = 0; i < ts.size(); i++) {
    if (taxa[i]) {
      taxa_list.insert(i);
    }
  }
}

Clade& Clade::operator=(const Clade& other) {
  taxa_list = other.taxa_list;
  taxa = other.taxa;
  ts = other.ts;
  return *this;
}

size_t Clade::size() const { return taxa_list.size(); }




string Clade::str() const {
  stringstream ss;
  vector<string> strings;

  for (Taxon i : taxa_list) {
    strings.push_back(ts[i]);
  }

  sort(strings.begin(), strings.end());
  
  ss << '{';
  for (string s : strings) {
    ss << s << ", ";
  }
  ss << '}';
  return ss.str();
}

string Clade::newick_str(TripartitionScorer& scorer, vector<Clade>& clades) {
  //  BOOST_LOG_TRIVIAL(info) << str() << endl;
  if (size() == 0) {
    return "";
  }
  if (size() == 1) {
    return ts[*taxa_list.begin()];
  }
  if (size() == 2) {
    stringstream ss;

    vector<Taxon> tv(taxa_list.begin(), taxa_list.end());
    
    ss << "("<<ts[tv[0]] << "," << ts[tv[1]] << ")" ;    
    
    return ss.str();

  }


  stringstream ss;
  pair<clade_bitset, clade_bitset>& subclades = scorer.get_subclades(taxa, clades);

  Clade c1(ts, subclades.first);
  Clade c2(ts, subclades.second);
  
  BOOST_LOG_TRIVIAL(debug) << str() << c1.str() << c2.str() << (int)scorer.get_score(taxa) <<endl;

  Tripartition tp(ts, *this, c1);
  
  ss << "(" << c1.newick_str(scorer, clades) << "," << c2.newick_str(scorer, clades) << ")";
  return ss.str();
}

void Clade::test() {
  TaxonSet ts;
  string str = string("{tx1, tx8, tx3, tx2, tx4}");
  cout << Clade(ts, str).str() << endl;
  cout << ts.str() << endl;
}

bool Clade::contains(const Clade& other) const {
  return (other.taxa & taxa) == other.taxa;
}
bool Clade::contains(const Taxon taxon) const {
  return taxa[taxon];
}

void Clade::add(const Taxon taxon) {
  taxa.set(taxon);
  taxa_list.insert(taxon);
}

Clade Clade::complement() const {
  Clade c(ts, ts.taxa_bs ^ taxa);
  return c;
}

Clade Clade::minus(const Clade& other) const {
  Clade c(ts, other.taxa ^ taxa);
  return c;
}


double Clade::score(TripartitionScorer& scorer, vector<Clade>& clades, unordered_set<clade_bitset>& cladetaxa) {
  double value;

  
  if (size() == 1) {
    value = 0;
    Clade eclade (ts);
    scorer.set_score(taxa, value, taxa, eclade.taxa);
    return value;
  }

  value = scorer.get_score(taxa);
  if (value != numeric_limits<double>::infinity()) {
    return value;
  }
  clade_bitset sub1, sub2;
  
  if (size() == 2) {
    Clade c1(ts);
    c1.add(*taxa_list.begin());
    Tripartition tp(ts, *this, c1);
    value = scorer.score(tp);
    scorer.set_score(taxa, value, tp.a1.taxa, tp.a2.taxa);
  }
  else {
    //    BOOST_LOG_TRIVIAL(debug) << "STARTING " << taxa_list.size() << endl << endl;
    for (Clade& subclade: clades) {
      if (!contains(subclade) || subclade.size() == 0 || subclade.size() >= size())
	continue;
      Tripartition tp(ts, *this, subclade);

      if (cladetaxa.count(tp.a1.taxa) == 0 || cladetaxa.count(tp.a2.taxa) == 0)
	continue;
      
      double score = scorer.score(tp) + tp.a1.score(scorer, clades, cladetaxa) + tp.a2.score(scorer, clades, cladetaxa);
      if (value == numeric_limits<double>::infinity() || (score < value) ) {
	value = score;
	sub1 = tp.a1.taxa;
	sub2 = tp.a2.taxa;	
      }
    }
    scorer.set_score(taxa, value, sub1, sub2);    
  }
  

  
  return value;
  
}

Tripartition::Tripartition(TaxonSet& ts, Clade& clade, Clade& subclade) :
  a1(ts), a2(ts), rest(ts)
{
  assert(clade.contains(subclade));
  a1 = clade.minus(subclade);
  a2 = subclade;
  rest = clade.complement();
}
