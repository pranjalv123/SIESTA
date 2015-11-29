#ifndef ASTRAL_INTERFACE_HPP__
#define ASTRAL_INTERFACE_HPP__

#include "Logger.hpp"

class AstralInterface {
private:
  string astralPath;
public:
  AstralInterface(string astralPath) : astralPath(astralPath) {}
  string getClades(string genetreefile) {
    string s = "java -jar " + astralPath + " -i " + genetreefile + " -k searchspace_norun -o /dev/null";
    if (!Logger::isEnabled("DEBUG") )
	s += " 2> /dev/null";

    PROGRESS << "Running ASTRAL to get clade set" << endl;
    DEBUG << "Using command line " << s << endl;
    
    FILE* stream = popen(s.c_str(), "r");
    char buffer[128];
    stringstream result;
    while(!feof(stream)) {
      if(fgets(buffer,128,stream) != NULL) {
	result << buffer;
      }
    }
    return result.str();
  }

  string getClades_exact(string genetreefile) {
    string s = "java -jar " + astralPath + " -x -i " + genetreefile + " -k searchspace_norun -o /dev/null";

    if (!Logger::isEnabled("DEBUG") )
	s += " 2> /dev/null";

    PROGRESS << "Running ASTRAL in exact mode to get clade set" << endl;
    DEBUG << "Using command line " << s << endl;
    
    FILE* stream = popen(s.c_str(), "r");
    char buffer[128];
    stringstream result;
    while(!feof(stream)) {
      if(fgets(buffer,128,stream) != NULL) {
	result << buffer;
      }
    }
    return result.str();
  }

  
  string getClades_limited(string genetreefile) {
    string s = "java -jar " + astralPath + " -p0 -i " + genetreefile + " -k searchspace_norun -o /dev/null";

    if (!Logger::isEnabled("DEBUG") )
	s += " 2> /dev/null";

    PROGRESS << "Running ASTRAL in limited mode to get clade set" << endl;
    DEBUG << "Using command line " << s << endl;
    
    FILE* stream = popen(s.c_str(), "r");
    char buffer[128];
    stringstream result;
    while(!feof(stream)) {
      if(fgets(buffer,128,stream) != NULL) {
	result << buffer;
      }
    }
    return result.str();
  }
};

#endif
