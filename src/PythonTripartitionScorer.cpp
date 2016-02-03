#include "TripartitionScorer.hpp"
#include "PythonTripartitionScorer.hpp"
#include "Options.hpp"
#include <python2.7/Python.h>

DEF_SCORER(PythonTripartitionScorer);

PythonTripartitionScorer::PythonTripartitionScorer(TaxonSet& ts) :
  TripartitionScorer(ts)
{
  string pfile;
  
  Py_Initialize();
  
  assert(Options::get("p pythonfile", &pfile));
  DEBUG << "Python module: " << pfile << endl;
  pName = PyString_FromString(pfile.c_str());
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  assert(pModule);

  pInitFn = PyObject_GetAttrString(pModule, "init");
  pScoreFn = PyObject_GetAttrString(pModule, "score");

  assert(pScoreFn);

  if (pInitFn && PyCallable_Check(pInitFn)) {
    DEBUG << "Executing init function" << endl;
    pArgs = PyTuple_New(1);
    
    PyTuple_SetItem(pArgs, 0, PyString_FromString(ts.str().c_str()));
    
    //    PyObject_CallObject(pInitFn, pArgs);
    PyObject_CallFunction(pInitFn, NULL);
    if (PyErr_Occurred())
       PyErr_Print();
    Py_DECREF(pArgs);
    Py_Finalize();
    exit(0);
  }
  else {
    ERR << "No init function" << endl;
  }
  Py_DECREF(pInitFn);
}

double PythonTripartitionScorer::score(const Tripartition& t) {
  return 0;
}
