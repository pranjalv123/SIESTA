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
  Options::get("p pythonfile", &pfile);
  DEBUG << "Python module: " << pfile << endl;
  pName = PyString_FromString(pfile.c_str());
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  if (PyErr_Occurred())
    PyErr_Print();
  
  assert(pModule);

  pInitFn = PyObject_GetAttrString(pModule, "init");
  pAdjustFn = PyObject_GetAttrString(pModule, "adjust");
  pScoreFn = PyObject_GetAttrString(pModule, "score");
  if (PyErr_Occurred())
    PyErr_Print();

  
  assert(pScoreFn);

  if (pInitFn && PyCallable_Check(pInitFn)) {
    DEBUG << "Executing init function" << endl;


    PyObject *pArgs, *pTaxonSet;

    PyObject *dendropy = PyImport_Import(PyString_FromString("dendropy"));
    assert(dendropy);

    PyObject *tn = PyObject_GetAttrString(dendropy, "TaxonNamespace");
    assert(tn);
    
    tn = PyObject_CallObject(tn, NULL);
    assert(tn);
    
    PyObject *add = PyObject_GetAttrString(tn, "new_taxon");
    assert(add);
    
    for (int i = 0; i < ts.size(); i++) {
      pArgs = PyTuple_New(1);
      PyTuple_SetItem(pArgs, 0, PyString_FromString(ts[i].c_str()));
      PyObject_CallObject(add, pArgs);
      Py_DECREF(pArgs);
    }
      

    pArgs = PyTuple_New(2);
    
    PyTuple_SetItem(pArgs, 0, tn);

    PyObject *ls = PyList_New(0);

    for (string& s : Options::argv) {
      PyList_Append(ls, PyString_FromString(s.c_str()));
    }
    
    PyTuple_SetItem(pArgs, 1, ls);
    
    PyObject_CallFunction(pInitFn, "O", pArgs);
    
    if (PyErr_Occurred())
       PyErr_Print();
    Py_DECREF(pArgs);
  }
  else {
    ERR << "No init function" << endl;
  }
  Py_DECREF(pInitFn);
}

double PythonTripartitionScorer::score(const Tripartition& t) 
{
  double output;
#pragma omp critical
  {
  PyObject *pArgs = PyTuple_New(3);
  PyObject *val1 = _PyLong_FromByteArray((unsigned char*)t.a1.taxa.data, t.a1.taxa.cap*sizeof(*t.a2.taxa.data), 1, 0);
  PyTuple_SetItem(pArgs, 0, val1);
  PyObject *val2 = _PyLong_FromByteArray((unsigned char*)t.a2.taxa.data, t.a2.taxa.cap*sizeof(*t.a2.taxa.data), 1, 0);
  PyTuple_SetItem(pArgs, 1, val2);
  PyObject *val3 = _PyLong_FromByteArray((unsigned char*)t.rest.taxa.data, t.rest.taxa.cap*sizeof(*t.a2.taxa.data), 1, 0);
  PyTuple_SetItem(pArgs, 2, val3);
  
  PyObject* retval = PyObject_CallFunction(pScoreFn, "O", pArgs);
  if (PyErr_Occurred()) {
    PyErr_Print();
    exit(1);
  }

  output = PyFloat_AsDouble(retval);
  
  Py_DECREF(val1);
  Py_DECREF(val2);
  Py_DECREF(val3);
  //Py_DECREF(pArgs);
  //Py_DECREF(retval);
  }
  return output;
    
}

double PythonTripartitionScorer::adjust_final_score(double d) {
  if (!pAdjustFn)
    return d;
  PyObject* pArgs = PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, PyFloat_FromDouble(d));
  PyObject* retval = PyObject_CallFunction(pAdjustFn, "O", pArgs);
  if (PyErr_Occurred()) {
    PyErr_Print();
    exit(1);
  }
  return PyFloat_AsDouble(retval);
  
}
