#include <sstream>
#include <cassert>
#include <ifstream>
#include "Quartet.hpp"

Quartet::Quartet(Taxon a, Taxon b, Taxon c, Taxon d)
{
  taxa[0] = a;
  taxa[1] = b;
  taxa[2] = c;
  taxa[3] = d;
}

double Quartet::parse(char* str){
  char* p = str;
  while(*p && *p != ':'){ p++; }
  assert(*p);
  *p = '\0';
  double weight = atof(p+1);
  
  if (str[0] == '(') {
    parse_newick(str);
  } else {
    parse_wqmc(str);
  }

  return weight;
}

void Quartet::parse_newick(char* c) {
  char* saveptr;
  int i = 0;
  while(token=strtok_r(c, "(),", &saveptr)) {
    cout << token << endl;
    c = NULL;
    taxa[i] = ts[token];
    i++;
  }
}


void Quartet::parse_wqmc(char* c) {
  char* saveptr;
  int i = 0;
  while(token=strtok_r(c, "|,", &saveptr)) {
    cout << token << endl;
    c = NULL;
    taxa[i] = ts[token];
    i++;
  }
}

string Quartet::str() {
  stringstream ss;
  ss << "((" << ts[a] << ", " << ts[b] << "),(" << ts[c] << ", " << ts[d] << "))";
  return ss.str();
}

QuartetDict::QuartetDict(TaxonSet& ts, string quartetfile) :
  ts(ts)
{
  Quartet q(ts);
  string s;
  double w;
  ifstream infile(quartetfile);

  while(!infile.eof()) {
    getline(infile, s);
    w = q.parse(s);
    array[q.taxa[0]][q.taxa[1]][q.taxa[2]][q.taxa[3]] = w;
    array[q.taxa[1]][q.taxa[0]][q.taxa[2]][q.taxa[3]] = w;
    array[q.taxa[0]][q.taxa[1]][q.taxa[3]][q.taxa[2]] = w;
    array[q.taxa[1]][q.taxa[0]][q.taxa[3]][q.taxa[2]] = w;
  }  
}

double operator[](Taxon a, Taxon b, Taxon c, Taxon d) {
  return array[a][b][c][d];
}


double operator[](Quartet& q) {
  return array[q.taxa[0]][q.taxa[1]][q.taxa[2]][q.taxa[3]];
}

string QuartetDict::str() {
  stringstream ss;
  for (int i = 0; i < ts.size(); i++) {
     for (int j = 0; j < i; j++) {
       for (int k = 0; k < j; k++) {
	 for (int l = 0; l < k; l++) {
	   Quartet q(i,j,k,l);
	   ss << q.str() << ":" << (*this)(q) << endl;
	 }
       }
     } 
  }
  return ss.str() ;
}

void QuartetDict::test() {
  TaxonSet ts;
  QuartetDict qd(ts, "quartetdict_test");
  cout << qd.str();
}


void Quartet::test() {
  TaxonSet ts;
  Quartet q(ts);
  q.parse("((1,2),(3,4))");
  assert(q.taxa[0] == ts["1"]);
  assert(q.taxa[1] == ts["2"]);
  assert(q.taxa[0] == ts["3"]);
  assert(q.taxa[1] == ts["4"]);
  cout << q.str();
}
