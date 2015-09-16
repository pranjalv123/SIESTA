#ifndef QUARTET_COUNTER_HPP__
#define QUARTET_COUNTER_HPP__

#include <vector>
#include <string>

#include "Clade.hpp"
#include "Quartet.hpp"

class QuartetCounter {
private:
  vector<vector<Clade> > trees;
public:
  QuartetCounter(TaxonSet& ts, vector<string> treelist);
  QuartetCounter(TaxonSet& ts, vector<vector<Clade> > treelist) : trees(treelist) {}
  
};

#endif
