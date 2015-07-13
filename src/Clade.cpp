#include "Clade.hpp"
#include "TripartitionScorer.hpp"
#include <sstream>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <cstring>

Clade::Clade(TaxonSet& tx_, string& str) :
  tx(tx_)
{
  int i = 1;

  char* cladestr = &(str[1]);
  char* token;
  char* saveptr;

  while(token = strtok_r(cladestr, ",} ", &saveptr)) {
    cladestr = NULL;
    add(tx[string(&(token[0]))]);
  }
}


Clade::Clade(TaxonSet& tx_) :
  tx(tx_)
{}

Clade::Clade(TaxonSet& tx_, bitset<128> taxa) :
  tx(tx_),
  taxa(taxa)
{
  for (int i = 0; i < tx.size(); i++) {
    if (taxa[i]) {
      taxa_list.insert(i);
    }
  }
}




string Clade::str() const {
  stringstream ss;
  vector<string> strings;

  for (Taxon i : taxa_list) {
    strings.push_back(tx[i]);
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
  cout << str() << endl;
  if (size() == 0) {
    return "";
  }
  if (size() == 1) {
    return tx[*taxa_list.begin()];
  }
  if (size() == 2) {
    stringstream ss;

    vector<Taxon> tv(taxa_list.begin(), taxa_list.end());
    
    ss << "("<<tx[tv[0]] << "," << tx[tv[1]] << ")" ;    
    
    return ss.str();

  }


  stringstream ss;
  pair<Clade, Clade>& subclades = scorer.get_subclades(taxa, clades);
  ss << "(" << subclades.first.newick_str(scorer, clades) << "," << subclades.second.newick_str(scorer, clades) << ")";
  return ss.str();
}

void Clade::test() {
  TaxonSet tx;
  string str = string("{tx1, tx8, tx3, tx2, tx4}");
  cout << Clade(tx, str).str() << endl;
  cout << tx.str() << endl;
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
  Clade c(tx, tx.taxa_bs ^ taxa);
  return c;
}

Clade Clade::minus(const Clade& other) const {
  Clade c(tx, other.taxa ^ taxa);
  return c;
}

int cladescorecount = 0;

double Clade::score(TripartitionScorer& scorer, vector<Clade>& clades, unordered_set<bitset<128> >& cladetaxa) {
  double value;

  
  if (size() == 1) {
    value = 0;
    return value;
  }

  value = scorer.get_score(taxa);
  if (value != numeric_limits<double>::infinity()) {
    return value;
  }

  
  if (size() == 2) {
    Clade c1(tx);
    c1.add(*taxa_list.begin());
    Tripartition tp(tx, *this, c1);
    value = scorer.score(tp);
    scorer.set_score(taxa, value, tp.a1, tp.a2);
  }

  else {
    
    cladescorecount++;
    if (cladescorecount % 1000 == 0) {
      //      cout << cladescorecount << "/" << clades.size() << endl;
    }
    for (Clade& subclade: clades) {
      if (!contains(subclade) || subclade.size() == 0 || subclade.size() >= size())
	continue;
      Tripartition tp(tx, *this, subclade);

      if (cladetaxa.count(tp.a1.taxa) == 0 || cladetaxa.count(tp.a2.taxa) == 0)
	continue;
      
      double score = scorer.score(tp) + tp.a1.score(scorer, clades, cladetaxa) + tp.a2.score(scorer, clades, cladetaxa);
      if (value == numeric_limits<double>::infinity() || score < value) {
	value = score;
	scorer.set_score(taxa, value, tp.a1, tp.a2);
      }
    }
    
  }
  
  return value;
  
}
