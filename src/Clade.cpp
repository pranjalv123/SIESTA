#include "Clade.hpp"

#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>

Clade::Clade(TaxonSet& tx_, string& str) :
  value(-1),
  best_subclade(NULL),
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
  value(-1),
  best_subclade(NULL),
  tx(tx_)
{
}




string Clade::str() {
  stringstream ss;
  ss << '{';
  for (int i = 0; i < taxa.size(); i++) {
    if (taxa[i]) 
      ss << tx[i] << ", ";
  }
  ss << '}';
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
  Clade c(tx);
  c.taxa = tx.taxa_bs ^ taxa;
  for (int i = 0; i < tx.size(); i++) {
    if (c.taxa[i]) {
      c.taxa_list.insert(i);
    }
  }
  return c;
}

Clade Clade::minus(const Clade& other) const {
  Clade c(tx);
  c.taxa = other.taxa ^ taxa;
  for (int i = 0; i < tx.size(); i++) {
    if (c.taxa[i]) {
      c.taxa_list.insert(i);
    }
  }
  return c;
}

double Clade::score(TripartitionScorer& scorer,  vector<Clade>& clades) {
  if (value != -1) {
    return value;
  }

  if (size() == 1) {
      clade.value = 0;
  }
  
  if (size() == 2) {
    Clade c1;
    Clade c2;
    c1.add(*taxa_list.begin());
    c2.add(*taxa_list.end());
    clade.value = scorer.score(Tripartition(tx, c1, c2)); 
  }

  else {
    for (Clade& subclade: clades) {
      if (!clade.includes(subclade))
	continue;
      
    }

  }

  return value;
  
}
