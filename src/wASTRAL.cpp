#include "Clade.hpp"
#include "Quartet.hpp"
#include "CladeSelector.hpp"
#include "TripartitionScorer.hpp"
#include "AstralInterface.hpp"
#include "Options.hpp"

#include <fstream>
#include <iostream>

#ifdef ENABLE_PROFILING
#include <gperftools/profiler.h>
#endif

void test() {
  Clade::test();
  Quartet::test();
  QuartetDict::test();
}

int main(int argc, char** argv) {
  TaxonSet ts;

  Options opt(argc, argv);
  
  
  if (opt.help || argc == 1) {
    cout << opt.desc << "\n";
    return 1;
  }

  if (opt.minimize && opt.maximize) {
    cerr << "ERROR: --minimize and --maximize are not compatible." << endl;
    return 1;
  }
  
  vector<Clade> clades;
  unordered_set<clade_bitset > cladetaxa;
  
  if (opt.cladefile.size()) {
    ifstream cladeFile(opt.cladefile);
    string s;
    while(!cladeFile.eof()) {
      getline(cladeFile, s);
      if (s.size() == 0)
	continue;
      clades.emplace_back(ts, s);
      cladetaxa.insert(clades.back().taxa);
    }    
  } 
  else if (opt.astralfile.size()) {
   AstralInterface ai(opt.astralfile);
    string cladesstr;
    if (opt.exact) {
      cladesstr = ai.getClades_exact(opt.genetreesfile, opt.verbose);
    } else {
      cladesstr = ai.getClades(opt.genetreesfile, opt.verbose);
    }
    stringstream cladess(cladesstr);
    string s;
    while(!cladess.eof()) {
      getline(cladess, s);
      if (s.size() == 0)
	continue;
      clades.emplace_back(ts, s);
      cladetaxa.insert(clades.back().taxa);
    }
  }

  string quartetFile(opt.quartetsfile);
  
  Clade alltaxa(ts);
  for (int i = 0; i < ts.size(); i++) {
    alltaxa.add(i);
  }

  clades.push_back(alltaxa);
  
  
  int count=0;

  QuartetDict qd(ts, quartetFile, opt.maximize);

  DPTripartitionScorer scorer(ts, qd);
  
  CladeSelector cs(ts, scorer, clades, cladetaxa);

  if (opt.profile) {
#ifdef ENABLE_PROFILING 
    ProfilerStart(opt.profilefile.c_str());
#else
    cerr << "wASTRAL must be compiled with ENABLE_PROFILING=ON for profiling to work!" << endl;
    return 1;
#endif
  }
  cs.run(opt.maximize);

#ifdef ENABLE_PROFILING
  if (opt.profile)
    ProfilerStop();
#endif
  if(opt.outputfile.size()) {
    ofstream outfile(opt.outputfile);
    outfile << cs.newick_tree << ';' << endl;
  }
  
  
}
