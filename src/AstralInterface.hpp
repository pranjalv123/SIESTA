#ifndef ASTRAL_INTERFACE_HPP__
#define ASTRAL_INTERFACE_HPP__

class AstralInterface {
private:
  string astralPath;
public:
  AstralInterface(string astralPath) : astralPath(astralPath) {}
  string getClades(string genetreefile, bool verbose) {
    string s = "java -jar " + astralPath + " -i " + genetreefile + " -k searchspace_norun -o /dev/null";
    if (!verbose) {
      s += " 2> /dev/null";
    }

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

  string getClades_exact(string genetreefile, bool verbose) {
    string s = "java -jar " + astralPath + " -x -i " + genetreefile + " -k searchspace_norun -o /dev/null";
    if (!verbose) {
      s += " 2> /dev/null";
    }
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
