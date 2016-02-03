#ifndef PYTHON_TRIPARTITION_SCORER_HPP__
#define PYTHON_TRIPARTITION_SCORER_HPP__

#include "TripartitionScorer.hpp"

#include <python2.7/Python.h>

class PythonTripartitionScorer : public TripartitionScorer{
public:
  DEC_SCORER(PythonTripartitionScorer);
  PythonTripartitionScorer(TaxonSet& ts);
  virtual double score(const Tripartition& t);


private:  
  PyObject *pName, *pModule, *pScoreFn, *pInitFn;
  PyObject *pArgs, *pValue;

  
};



#endif // PYTHON_TRIPARTITION_SCORER_HPP__
