#include "wASTRAL.hpp"
#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv) {
  vector<string> wastral_args;

  string input;
  string output;
  string extra;
  int debug = 0;
  
  for(int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-i") {
      assert(argc > i+1);
      i++;
      input = string(realpath(argv[i], NULL));
    }
    if (string(argv[i]) == "-o") {
      assert(argc > i+1);
      i++;
      output = string(realpath(argv[i], NULL));
    }
    if (string(argv[i]) == "-e") {
      assert(argc > i+1);
      i++;
      extra = string(realpath(argv[i], NULL));
    }
    if (string(argv[i]) == "--debug") {
      debug=1;
    }
  }

  if (input.size() == 0) {
    cerr << "FastRFS -i <source tree file> [-o <output file>] [-e <extra trees>]" << endl;
    exit(-1);
  }
  
  wastral_args.push_back("-a");
  wastral_args.push_back("astral.4.7.8.jar");
  wastral_args.push_back("-c");
  wastral_args.push_back("FastRF");
  wastral_args.push_back("--maximize");
  
  wastral_args.push_back("-g");
  wastral_args.push_back(input);

  if (debug) {
    wastral_args.push_back("--verbose");
    wastral_args.push_back("debug");
  }
  
  if (extra.size()) {
    wastral_args.push_back("--extraextra");
    wastral_args.push_back("-e");
    wastral_args.push_back(extra);
  }
  
  if (output.size()) {
    wastral_args.push_back("-o");
    wastral_args.push_back(output);
  }

  vector<const char*> wastral_args_c;
  wastral_args_c.push_back(argv[0]);
  for (string& s : wastral_args) {
    wastral_args_c.push_back(s.c_str());
  }

  if (debug) {
    for (const char* s : wastral_args_c) {
      cout << s << " ";
    }
    cout << endl;
  }
  
  wASTRAL(wastral_args_c.size() , wastral_args_c.data());
}
