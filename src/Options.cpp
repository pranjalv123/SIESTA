#include "Options.hpp"
#include <getopt.h>

string Options::desc =
  "Options for wASTRAL:\n\
  -h [ --help ]          produce help message\n\
  -q [ --quartets ] arg  Input quartet weight file\n\
  -c [ --cladefile ] arg Input clades file\n\
  -g [ --genetrees ] arg Input gene trees, will generate clades with ASTRAL.\n\
                          The ASTRAL executable must be present in the current \n\
                         directory as astral.jar, or it can be specified with \n\
                         --astral/-a\n\
  -a [ --astral ] arg    Path to ASTRAL .jar file\n\
  -x [ --exact ] arg     Get exact solution\n\
  -o [ --output ] arg    Output tree file\n\
  --minimize             Minimize sum of quartet weights (default)\n\
  --maximize             Maximize sum of quartet weights\n\
  -s [ --score ] arg     Output the score of this tree\n\
  -v [ --verbose ]       Verbose output (debugging)\n\
  --profile              Enable profiling with pprof (debugging)\n\
" ;

Options::Options(int argc, char** argv) {
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"test", no_argument, &test, 1},

    
    {"astral", required_argument, 0, 'a'},
    {"exact", no_argument, 0, 'x'},
    
    {"quartets", required_argument, 0, 'q'},    
    {"cladefile", required_argument, 0, 'c'},
    {"genetrees", required_argument, 0, 'g'},
    {"score", required_argument, 0, 's'},
    {"output", required_argument, 0, 'o'},
    
    {"maximize", no_argument, &maximize, 1},
    {"minimize", no_argument, &minimize, 1},
    
    {"profile", required_argument, 0, 'p'}

  };

  string optstring(":hva:xq:c:g:s:o:p:");

  int opt = 0;
  int long_index = 0;
  verbose = false;
  exact = false;
  help = false;
  profile = false;
  maximize = false;
  minimize = false;
  test = false;
  while ((opt = getopt_long(argc, argv, optstring.c_str(), long_options, &long_index)) != -1) {
    switch(opt) {
    case 'v': verbose = true;
      break;
    case 'q': quartetsfile = optarg;
      break;
    case 'c': cladefile = optarg;
      break;
    case 'g': genetreesfile = optarg;
      break;
    case 'a': astralfile = optarg;
      break;
    case 'x': exact = true;
      break;
    case 'o': outputfile = optarg;
      break;
    case 'h': help = true;
      break;
    case 'p': profilefile = optarg; profile = true;
      break;
    }
  }
}
