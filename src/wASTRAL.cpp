#include "Clade.hpp"
#include "CladeExtractor.hpp"
#include "Quartet.hpp"
#include "CladeSelector.hpp"
#include "TripartitionScorer.hpp"
#include "AstralInterface.hpp"
#include "Options.hpp"
#include "Logger.hpp"

#include "RFTripartitionScorer.hpp"
#include "DPTripartitionScorer.hpp"
#include "BryantSteelTripartitionScorer.hpp"
#include "PythonTripartitionScorer.hpp"

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

  //parse the command-line arguments, set up the logger, set up profiling, ...
  
  Options::init(argc, argv);
  
  if(Options::get("h help") || argc==1 ){
    cout << " --a astral &astralfile \n\
--c criterion &heuristic\n\
--g genetrees &genetreesfile\n\
--g genetrees &treesfile\n\
--h help\n\
--maximize\n\
--maximize\n\
--o output &output\n\
--p pythonfile\n\
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
  
  bool profile = Options::get("profile", &profilefile);
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


  //This also initializes the clades and gets them from whatever source is specified
  //This could be gene trees with ASTRAL or clades that are explicitly given
  
  TaxonSet& ts = CladeExtractor::get_taxonset();


  //Create whatever tripartition scorer we want to use, e.g. RF, BryantSteel, etc.
  //See TripartitionScorer.hpp for these definitions
  TripartitionScorer* tps = TripartitionScorerFactory::createInstance(heuristic, ts);
  if (!tps)
    tps = TripartitionScorerFactory::createInstance(heuristic + "TripartitionScorer", ts);
  if (!tps) {
    ERR << "Invalid criterion: make sure you have supplied the -c argument\n";
    exit(-1);
      
  }

  vector<Clade> cladev(CladeExtractor::get_clades().begin(), CladeExtractor::get_clades().end());

  //Actually run the DP algorithm
  CladeSelector cs(ts, *tps, cladev, CladeExtractor::get_cladetaxa());

  bool maximize = Options::get("maximize");
  
  double score = cs.run(maximize);

  
#ifdef ENABLE_PROFILING
  if (profile)
    ProfilerStop();
#endif

  //Output results
  string output;
  if(Options::get("o output", &output)) {
    ofstream outfile(output);
    if (Options::get("s score")) {
      outfile << tps->adjust_final_score(score) << endl;
    } else {
      outfile << cs.newick_tree << ';' << endl;
    }
  }
  
  
}
