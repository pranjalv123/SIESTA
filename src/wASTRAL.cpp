#include "Clade.hpp"
#include "CladeExtractor.hpp"
#include "Quartet.hpp"
#include "CladeSelector.hpp"
#include "TripartitionScorer.hpp"
#include "AstralInterface.hpp"
#include "Options.hpp"
#include "Logger.hpp"
#include "wASTRAL.hpp"

#include "RFTripartitionScorer.hpp"
#include "FastRFTripartitionScorer.hpp"
#include "BryantSteelTripartitionScorer.hpp"
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

  double count=0, ntrees=0, defective=0;
  
  if (scoremat) {

    if (Logger::isEnabled("DEBUG") )
	for (size_t i = 0; i < cladev.size(); i++) {
	  for (size_t j = 0; j < cladev.size(); j++) {
	    cout << (*scoremat)[i][j] << "\t";	
	  }
	  cout << endl;
	}


   
    unordered_map<clade_bitset, double > count_cache;

    string defective_str;


    unordered_map< clade_bitset, unordered_map<double, double> > defective_cache;

    
    if (!Options::get("defective", &defective_str)){
       count = cladev.back().optimal_subtree_count(*tps, count_cache);
    }
    else {
      defective = atof(&(defective_str[0]));
      
      unordered_map<clade_bitset, int> clade_indices;    
  
      for ( size_t i = 0; i < cladev.size(); i++ ) {
	cladev[i].myIndex = i;
	clade_indices[cladev[i].taxa] = i;
      }
      for (int i = 0; i <= defective; i++) {
	count += cladev.back().defective_subtree_count(*tps, *scoremat, cladev, clade_indices, i, defective_cache);
      }

      
      for (Clade& c : cladev) {
	for (int i = 0; i <= defective; i++) {
	  count_cache[c.taxa] = 0;
	  DEBUG << c.str() << "\t" << i << "\t" << defective_cache[c.taxa][i] << endl;
	  count_cache[c.taxa] += defective_cache[c.taxa][i];
	}		
      }    
    }
    

    ntrees = count/(2 * cladev.back().size() - 3);    
    INFO << "Found " << ntrees << " optimal trees" << endl;
    sort(cladev.begin(), cladev.end(), [](const ScorableClade& a, const ScorableClade& b){ return a.size() > b.size(); });

    
    if (Options::get("listcladecounts")) {

      for (auto& c : cladev) {
	if (defective)
	  cout << c.str() << "\t" << c.appearances_in_defective_trees(*tps, defective, defective_cache) << "\t" << count_cache[c.taxa] << "\t" <<count_cache[c.complement().taxa] << endl;
	else
	  cout << c.str() << "\t" << c.appearances_in_optimal_trees(*tps, count_cache) << "\t" << count_cache[c.taxa] << "\t" <<count_cache[c.complement().taxa] << endl;
      }
    }

    string consensus_str;

    
    if (Options::get("consensus", &consensus_str)) {
      double consensus = atof(&(consensus_str[0])) * ntrees;
      
      vector<ScorableClade> consensusclades;
      consensusclades.push_back(cladev[0]);
      for (auto& c : cladev) {
	if ((defective && c.appearances_in_defective_trees(*tps, defective, defective_cache) >= consensus) || (c.appearances_in_optimal_trees(*tps, count_cache) >= consensus)) {
	  consensusclades.push_back(c);
	}
      }

      DEBUG << "Consensusclades" << consensusclades.size() << endl;
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
      outfile << ntrees << endl;
    }else {
      outfile << cs.newick_tree << ';' << endl;
    }
  }
  
  return 0;
}
