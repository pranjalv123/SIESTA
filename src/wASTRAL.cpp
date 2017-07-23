#include "Clade.hpp"
#include "CladeExtractor.hpp"
#include "Quartet.hpp"
#include "CladeSelector.hpp"
#include "TripartitionScorer.hpp"
#include "AstralInterface.hpp"
#include <util/Options.hpp>
#include <util/Logger.hpp>
#include "wASTRAL.hpp"

#include "RFTripartitionScorer.hpp"
#include "FastRFTripartitionScorer.hpp"
#include "BryantSteelTripartitionScorer.hpp"
#include "AstralTripartitionScorer.hpp"
#include "DPTripartitionScorer.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

#ifdef ENABLE_PROFILING
#include <gperftools/profiler.h>
#endif

void test() {
  CladeExtractor::test();
  Clade::test();
  Quartet::test();
  QuartetDict::test();
}



int wASTRAL(int argc, const char** argv) {

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
    ERR << TripartitionScorerFactory::instanceList();
    exit(-1);
      
  }

  vector<ScorableClade> cladev(CladeExtractor::get_clades().begin(), CladeExtractor::get_clades().end());


  sort(cladev.begin(), cladev.end(), [](const ScorableClade& a, const ScorableClade& b){ return a.size() < b.size(); });

  twod_mat* scoremat = 0;
  if (Options::get("matrix")) {
      scoremat = new twod_mat(boost::extents[cladev.size()][cladev.size()]);
      for (size_t i = 0; i < cladev.size(); i++) {
	for (size_t j = 0; j < cladev.size(); j++) {
	  (*scoremat)[i][j] = (double)nan("");
	}
      }
  }
  
  //Actually run the DP algorithm
  BasicCladeSelector cs(ts, *tps, cladev, CladeExtractor::get_cladetaxa());
  
  bool maximize = Options::get("maximize");
  
  double score = cs.run(maximize, scoremat);
  string consensus_tree;

  __int128 count=0, ntrees=0, defective=0;
  
  if (scoremat) {

    if (Logger::isEnabled("DEBUG") )
	for (size_t i = 0; i < cladev.size(); i++) {
	  for (size_t j = 0; j < cladev.size(); j++) {
	    cout << (*scoremat)[i][j] << "\t";	
	  }
	  cout << endl;
	}

  }
   
  unordered_map<clade_bitset, __int128 > count_cache;

  string defective_str;


  //  unordered_map< clade_bitset, unordered_map<double, double> > defective_cache;
    
    
  
  count = cladev.back().optimal_subtree_count(*tps, count_cache);
  
  ntrees = count/(2 * cladev.back().size() - 3);    
  cerr << "Found " << (double)ntrees << " optimal trees" << endl;
  sort(cladev.begin(), cladev.end(), [](const ScorableClade& a, const ScorableClade& b){ return a.size() > b.size(); });

    
  if (Options::get("listcladecounts")) {
    
    for (auto& c : cladev) {
      // if (defective)
      // 	cout << c.str() << "\t" << c.appearances_in_defective_trees(*tps, defective, defective_cache) << "\t" << count_cache[c.taxa] << "\t" <<count_cache[c.complement().taxa] << endl;
      // else
      //      cout << c.str() << "\t" <<(c.appearances_in_optimal_trees(*tps, count_cache) << "\t" << count_cache[c.taxa] << "\t" <<count_cache[c.complement().taxa] << endl;
    }
  }
  
  string consensus_str;
  
  
  if (Options::get("consensus", &consensus_str)) {
    __int128 consensus = atof(&(consensus_str[0])) * ntrees;    
    vector<ScorableClade> consensusclades;
    consensusclades.push_back(cladev[0]);
    
    sort(cladev.begin(), cladev.end(), [ tps, &count_cache ](const ScorableClade& a, const ScorableClade& b){ return a.appearances_in_optimal_trees(*tps, count_cache) > b.appearances_in_optimal_trees(*tps, count_cache); });
    

    

    for (auto& c : cladev) {
      //      cout << "\t" << c.appearances_in_optimal_trees(*tps, count_cache) << "\t" << consensus << "\t" << c.appearances_in_optimal_trees(*tps, count_cache)  - consensus << endl << "\t" << count_cache[c.taxa] - consensus << "\t" << count_cache[c.complement().taxa] - consensus << endl;
      if ((c.appearances_in_optimal_trees(*tps, count_cache) >= consensus)) {
	bool compat = true;
	for (auto& cc : consensusclades) {
	  //DEBUG << c.str() << "\t" << cc.str() << endl;
	  //DEBUG << c.taxa.cap << "\t" << cc.taxa.cap << endl;
	  int os = c.overlap_size(cc);
	  
	  if (!(os == 0 || os == c.size() || os == cc.size())) {
	    compat=false;
	    break;
	  }
	}
	if (compat)
	  consensusclades.push_back(c);
      }
    }

    sort(cladev.begin(), cladev.end(), [](const ScorableClade& a, const ScorableClade& b){ return a.size() > b.size(); });
    sort(consensusclades.begin(), consensusclades.end(), [](const ScorableClade& a, const ScorableClade& b){ return a.size() > b.size(); });
    
    INFO << "Found " << consensusclades.size() << " consensus clades" << endl;
    if (Logger::isEnabled("DEBUG") ) {
      for (auto& c : consensusclades)
	DEBUG << c.str() << endl;
    }
    
    vector<double> supports;
    supports.push_back(1.0);
    for (int i = 1; i < consensusclades.size(); i++) {
      supports.push_back(consensusclades[i].appearances_in_optimal_trees(*tps, count_cache)/ntrees);
    }
    
    consensus_tree = consensusclades[0].newick_str(consensusclades, supports, 0);
    cout << "Consensus tree: " << consensus_tree << endl;
    
      
  }
  
  
  if (Options::get("listalltrees")) {
    unordered_map<clade_bitset, vector<string> > strcache;
    for (ScorableClade& c : cladev) {
      c.all_newick_strs(*tps, cladev, *scoremat, strcache);
    }
    
    for (string s : cladev[0].all_newick_strs(*tps, cladev, *scoremat, strcache)) {
      cout << s << endl;      
    }
    
  }
    


  
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
    } else if (Options::get("consensus")) {
      outfile << consensus_tree << ';' << endl;
    } else if (Options::get("counttrees")) {
      outfile << (double)ntrees << endl;
    } else if (Options::get("listalltrees")) {
      unordered_map<clade_bitset, vector<string> > strcache;
      for (ScorableClade& c : cladev) {
	c.all_newick_strs(*tps, cladev, *scoremat, strcache);
      }
      
      for (string s : cladev[0].all_newick_strs(*tps, cladev, *scoremat, strcache)) {
	outfile << s << ";" << endl;      
      }
      
    } else {
      outfile << cs.newick_tree << ';' << endl;
    }
  }
  
  return 0;
}
