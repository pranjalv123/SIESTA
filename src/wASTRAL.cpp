#include "Clade.hpp"
#include "CladeExtractor.hpp"
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
  CladeExtractor::test();
  Clade::test();
  Quartet::test();
  QuartetDict::test();
}

int main(int argc, char** argv) {
  Options opt(argc, argv);


  if (opt.profile) {
#ifdef ENABLE_PROFILING 
    ProfilerStart(opt.profilefile.c_str());
#else
    cerr << "wASTRAL must be compiled with ENABLE_PROFILING=ON for profiling to work!" << endl;
    return 1;
#endif
  }
  
  
  if (opt.help || argc == 1) {
    cout << opt.desc << "\n";
    return 1;
  }

  if (opt.test) {
    cout << "TESTING" << endl;
    test();
    return 0;
  }
  
  if (opt.minimize && opt.maximize) {
    cerr << "ERROR: --minimize and --maximize are not compatible." << endl;
    return 1;
  }



  
  unordered_set<Clade> clade_set;

  unordered_set<clade_bitset > cladetaxa;

  stringstream clade_stream;
  
  if (opt.cladefile.size()) {
    ifstream cladeFile(opt.cladefile);
    clade_stream << cladeFile;
  }
  if (opt.astralfile.size()) {
    AstralInterface ai(opt.astralfile);	
    if (opt.exact) {
      clade_stream << ai.getClades_exact(opt.genetreesfile, opt.verbose);
    } else {
      clade_stream << ai.getClades(opt.genetreesfile, opt.verbose);
    }
  }

  TaxonSet ts(clade_stream.str());

  stringstream ss(clade_stream.str());
  string s;
  while (!ss.eof()) {
    getline(ss, s);
    Clade c(ts, s);
    clade_set.insert(c);
    cladetaxa.insert(c.taxa);
  }
  
  string quartetFile(opt.quartetsfile);
  
  Clade alltaxa(ts);
  for (size_t i = 0; i < ts.size(); i++) {
    alltaxa.add(i);
  }

  clade_set.insert(alltaxa);
  
  vector<Clade> clades(clade_set.begin(), clade_set.end());

  QuartetDict qd(ts, quartetFile, opt.maximize);

  //DPTripartitionScorer scorer(ts, qd);
  //BryantSteelTripartitionScorer scorer(ts, qd, clades);

  RFTripartitionScorer scorer(ts, opt.genetreesfile);

  
  
  CladeSelector cs(ts, scorer, clades, cladetaxa);

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
