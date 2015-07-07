#include "Clade.hpp"

#include <sstream>
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

    cout << token << endl;
    cladestr = NULL;
    taxa.push_back(tx[string(&(token[0]))]);
  }
  sort(taxa.begin(), taxa.end());
}

string Clade::str() {
  stringstream ss;
  ss << '{';
  for (int i : taxa) {
    ss << i << ", ";
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
