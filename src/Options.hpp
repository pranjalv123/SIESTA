#ifndef OPTIONS_HPP__
#define OPTIONS_HPP__
#include <string>
using namespace std;

class Options {
public:
  int verbose;
  string quartetsfile;
  string cladefile;
  string genetreesfile;
  string astralfile;
  int exact;
  string outputfile;
  int maximize;
  int minimize;
  string scoretree;
  bool profile;
  string profilefile;
  int help;
  int test;
  
  Options(int argc, char** argv);
  
  static string desc;
  
};

#endif
