#include "Clade.hpp"
#include "CladeExtractor.hpp"
#include "Quartet.hpp"
#include "CladeSelector.hpp"
#include "TripartitionScorer.hpp"
#include "AstralInterface.hpp"
#include "Options.hpp"
#include "Logger.hpp"

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
  Options::init(argc, argv);
  
  if(Options::get("h help") || argc==1 ){
    cout << " --a astral &astralfile \n\
--c criterion &heuristic\n\
--g genetrees &genetreesfile\n\
--g genetrees &treesfile\n\
--h help) || argc==1 \n\
--maximize\n\
--maximize\n\
--o output &output\n\
--profile\n\
--q quartets &quartetFile\n\
--s score &scoretree\n\
string opts string* arg\n\
--v verbose &level\n\
--X cladefile &cladefile\n\
--x exact" << endl;
    exit(0);
  }

  Logger::get();
  
  string profilefile;
  
  bool profile = Options::get("profile");
  string heuristic = "DPTripartitionScorer";
  Options::get("c criterion", &heuristic);
  
  INFO << "Using heuristic " << heuristic << endl;
  
  
  if (profile) {
#ifdef ENABLE_PROFILING 
    ProfilerStart(profilefile.c_str());
#else
    cerr << "wASTRAL must be compiled with ENABLE_PROFILING=ON for profiling to work!" << endl;
    return 1;
#endif
  }  

  TaxonSet& ts = CladeExtractor::get_taxonset();
  
  TripartitionScorer* tps = TripartitionScorerFactory::createInstance(heuristic, ts);

  vector<Clade> cladev(CladeExtractor::get_clades().begin(), CladeExtractor::get_clades().end());
  
  CladeSelector cs(ts, *tps, cladev, CladeExtractor::get_cladetaxa());

  bool maximize = Options::get("maximize");
  
  cs.run(maximize);

  
#ifdef ENABLE_PROFILING
  if (profile)
    ProfilerStop();
#endif

  string output;
  if(Options::get("o output", &output)) {
    ofstream outfile(output);
    outfile << cs.newick_tree << ';' << endl;
  }
  
  
}
