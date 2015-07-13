#include "Clade.hpp"
#include "Quartet.hpp"
#include "CladeSelector.hpp"
#include "TripartitionScorer.hpp"

#include <fstream>
#include <google/profiler.h>


void test() {
  Clade::test();
  Quartet::test();
  QuartetDict::test();
}

int main(int argc, char** argv) {
  TaxonSet ts;
  
  ifstream cladeFile(argv[1]);
  string quartetFile(argv[2]);

  vector<Clade> clades;
  unordered_set<bitset<128> > cladetaxa;
  string s;
  while(!cladeFile.eof()) {
    getline(cladeFile, s);
    if (s.size() == 0)
      continue;
    clades.emplace_back(ts, s);
    cladetaxa.insert(clades.back().taxa);
  }
  Clade alltaxa(ts);
  for (int i = 0; i < ts.size(); i++) {
    alltaxa.add(i);
  }

  clades.push_back(alltaxa);
  
  cout << ts.str() << endl;
  
  int count=0;
  
  QuartetDict qd(ts, quartetFile);

  DPTripartitionScorer scorer(ts, qd);
  
  CladeSelector cs(ts, scorer, clades, cladetaxa);
  //  ProfilerStart("prof");
  cs.run();
  //  ProfilerStop();
}
