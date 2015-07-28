#include "Clade.hpp"
#include "Quartet.hpp"
#include "CladeSelector.hpp"
#include "TripartitionScorer.hpp"
#include "AstralInterface.hpp"

#include <fstream>
#include <google/profiler.h>
#include <boost/program_options.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace po = boost::program_options;

void test() {
  Clade::test();
  Quartet::test();
  QuartetDict::test();
}

int main(int argc, char** argv) {
  TaxonSet ts;
  
  po::options_description desc("Options for wASTRAL");
  desc.add_options()
    ("help,h", "produce help message")
    ("quartets,q", po::value<string>(), "Input quartet weight file")
    ("cladefile,c", po::value<string>(), "Input clades file")
    ("genetrees,g", po::value<string>(), "Input gene trees, will generate clades with ASTRAL.\n The ASTRAL executable must be present in the current directory as astral.jar, or it can be specified with --astral/-a")
    ("astral,a", po::value<string>(), "Path to ASTRAL .jar file")
    ("exact,x", po::value<string>(), "Get exact solution")
    ("output,o", po::value<string>(), "Output tree file")
    ("minimize", "Minimize sum of quartet weights (default)")
    ("maximize", "Maximize sum of quartet weights")
    ("score,s", po::value<string>(), "Output the score of this tree")
    ("verbose,v", "Verbose output (debugging)")
    ("profile", "Enable profiling with pprof (debugging)")
    ;
  

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help") || argc == 1) {
    BOOST_LOG_TRIVIAL(warning) << desc << "\n";
    return 1;
  }

  if (!vm.count("verbose")) {
    boost::log::core::get()->set_filter
    (
     boost::log::trivial::severity >= boost::log::trivial::info
    );
  }

  if (vm.count("minimize") && vm.count("maximize")) {
    BOOST_LOG_TRIVIAL(error) << "ERROR: --minimize and --maximize are not compatible." << endl;
    return 1;
  }

  vector<Clade> clades;
  unordered_set<clade_bitset > cladetaxa;
  
  if (vm.count("cladefile")) {
    ifstream cladeFile(vm["cladefile"].as<string>());
    string s;
    while(!cladeFile.eof()) {
      getline(cladeFile, s);
      if (s.size() == 0)
	continue;
      clades.emplace_back(ts, s);
      cladetaxa.insert(clades.back().taxa);
    }    
  } 
  else if (vm.count("astral")) {
    AstralInterface ai(vm["astral"].as<string>());
    string cladesstr;
    if (vm.count("exact")) {
      cladesstr = ai.getClades_exact(vm["genetrees"].as<string>(), vm.count("verbose"));
    } else {
      cladesstr = ai.getClades(vm["genetrees"].as<string>(), vm.count("verbose"));
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
  
  string quartetFile(vm["quartets"].as<string>());
  
  Clade alltaxa(ts);
  for (int i = 0; i < ts.size(); i++) {
    alltaxa.add(i);
  }

  clades.push_back(alltaxa);
  
  BOOST_LOG_TRIVIAL(debug) << ts.str() << endl;
  
  int count=0;
  
  QuartetDict qd(ts, quartetFile, vm.count("maximize"));

  DPTripartitionScorer scorer(ts, qd);
  
  CladeSelector cs(ts, scorer, clades, cladetaxa);
  if (vm.count("profile"))
    ProfilerStart(vm["profile"].as<string>().c_str());

  cs.run(vm.count("maximize"));

  if (vm.count("profile"))
      ProfilerStop();

  if(vm.count("output")) {
    ofstream outfile(vm["output"].as<string>());
    outfile << cs.newick_tree << ';' << endl;
  }

  
}
