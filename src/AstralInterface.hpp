#ifndef ASTRAL_INTERFACE_HPP__
#define ASTRAL_INTERFACE_HPP__

#include <util/Logger.hpp>
#include <newick.hpp>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>

class AstralInterface {
private:
  string astralPath;
public:
   string remapped_treefile(string& input, TaxonSet& ts) {
    string line;
    ifstream ifile(input);
    
    char name[] = "/tmp/fileXXXXXX";
    mkstemp(name);

    ofstream of(name);
    while(!ifile.eof()) {
      getline(ifile, line);
      if (line.size())
	of << map_newick_names(line, ts) << endl;
    }
    
    return name;
  }
  
  AstralInterface(string astralPath) : astralPath(astralPath) {}
  string getClades(string genetreefile, string extratreesfile, bool exact, bool limited) {
    unordered_set<string> taxa;
    char buffer[128];
    stringstream gtrees_stream;


    FILE* gtrees_file_stream = fopen(genetreefile.c_str(), "r");

    assert(gtrees_file_stream);
    
    while (!feof(gtrees_file_stream)) {
      if (fgets(buffer, 128, gtrees_file_stream) != NULL) {
	gtrees_stream << buffer;
      }
    }

    int ntaxa = newick_to_ts(gtrees_stream.str(), taxa);

    
    TaxonSet ts(ntaxa);
    
    string s = "java -jar " + astralPath + " -i " + remapped_treefile(genetreefile, ts) + " -k searchspace_norun -o /dev/null";
    if (exact) {
      s += " -x ";
    }
    if (limited) {
      s += " -p0 ";
    }
    if (extratreesfile.size()) {
      s += " -e " + remapped_treefile(extratreesfile, ts);
    }
    if (!Logger::isEnabled("DEBUG") )
	s += " 2> /dev/null";

    PROGRESS << "Running ASTRAL to get clade set" << endl;
    PROGRESS << "Using command line " << s << endl;

    
    
    FILE* stream = popen(s.c_str(), "r");

    stringstream result;
    while(!feof(stream)) {
      if(fgets(buffer,128,stream) != NULL) {
	result << buffer;
      }
    }

    stringstream cladestream_mapped(result.str());
    string line;

    stringstream unmapped;
    
    while (!cladestream_mapped.eof()) {
      getline(cladestream_mapped, line);
      unmapped << unmap_clade_names(line, ts) << endl;
    }

    
    return unmapped.str();
  }

};

#endif
