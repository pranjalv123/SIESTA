#include "Clade.hpp"
#include "Options.hpp"
#include "Logger.hpp"
#include "TripartitionScorer.hpp"
#include <sstream>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>

Clade::Clade(TaxonSet& ts_, string& str) :
  taxa(ts_.size()),
  ts(ts_)
{
  char* cladestr = &(str[1]);
  char* token;
  char* saveptr;

  while((token = strtok_r(cladestr, ",} ", &saveptr))) {
    cladestr = NULL;
    add(ts[string(&(token[0]))]);
  }
}

Clade::Clade(TaxonSet& ts_) :
  taxa(ts_.size()),
  ts(ts_)
{}

Clade::Clade(TaxonSet& ts_, clade_bitset& taxa) :
  taxa(taxa),
  ts(ts_)
{
}

Clade::Clade(TaxonSet& ts_, unordered_set<Taxon>& taxa) :
  taxa(ts_.size()),
  ts(ts_)
{
  for (Taxon t : taxa) {
    add(t);
  }
}

Clade::Clade(const Clade& other) :
  taxa(other.taxa),
  ts(other.ts)
{
}



Clade& Clade::operator=(const Clade& other) {
   taxa = other.taxa;
   ts = other.ts;
   return *this;
}

bool Clade::operator==(const Clade& other) const {
   return taxa == other.taxa;
}





string Clade::str() const {
  stringstream ss;
  vector<string> strings;

  for (Taxon i : *this) {
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

  Clade c1(ts, subclades.first);
  Clade c2(ts, subclades.second);
  
  //  BOOST_LOG_TRIVIAL(debug) << str() << c1.str() << c2.str() << (int)scorer.get_score(taxa) <<endl;

  Tripartition tp(ts, *this, c1);
  
  ss << "(" << c1.newick_str(scorer, clades) << "," << c2.newick_str(scorer, clades) << ")";
  return ss.str();
}

void Clade::test() {
  string str = string("{tx1, tx8, tx3, tx2, tx4}");
  TaxonSet ts(str);
  cout << Clade(ts, str).str() << endl;
  cout << ts.str() << endl;
}

bool Clade::contains(const Clade& other) const {
  BitVectorFixed overlap(other.taxa & taxa);
  return overlap == other.taxa;
}
bool Clade::contains(const Taxon taxon) const {
  return taxa.get(taxon);
}

void Clade::add(const Taxon taxon) {
  taxa.set(taxon);
}

Clade Clade::complement() const {
  BitVectorFixed comp = ts.taxa_bs & (~taxa);
  Clade c(ts, comp);
  return c;
}

Clade Clade::minus(const Clade& other) const {
  BitVectorFixed m(taxa & (~other.taxa));
  Clade c(ts, m);
  return c;
}

int Clade::size() const {
  return taxa.popcount();
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
  if (!std::isnan(value)) {
    return value;
  }
  clade_bitset sub1(ts.size()), sub2(ts.size());

  int invert = 1;
  if (Options::get("maximize")) {
    invert = -1;
  }
  
  if (size() == 2) {
    Clade c1(ts);
    c1.add(*begin());
    Tripartition tp(ts, *this, c1);
    value = invert * scorer.score(tp);

    scorer.set_score(taxa, value, tp.a1.taxa, tp.a2.taxa);
  }
  else {
    for (Clade& subclade: clades) {
      if (subclade.size() >= size() || !contains(subclade) || subclade.size() == 0 )
	continue;

      Tripartition tp(ts, *this, subclade);

      if (cladetaxa.count(tp.a1.taxa) == 0 || cladetaxa.count(tp.a2.taxa) == 0)
	continue;

      
      double score = invert * scorer.score(tp) + tp.a1.score(scorer, clades, cladetaxa) + tp.a2.score(scorer, clades, cladetaxa);
      if (std::isnan(value) || (score < value) ) {
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
  a1(clade.minus(subclade)),
  a2(subclade),
  rest(clade.complement())
{
  assert(clade.contains(subclade));
}

void Clade::do_swap(Clade& other) {
  std::swap(taxa, other.taxa);
}


namespace std
{
    template<>
    void swap<Clade>(Clade& lhs, Clade& rhs)
    {
      lhs.do_swap(rhs);
    }
}
