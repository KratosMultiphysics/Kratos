#include <Python.h>
#include <iostream>

int solve(/*std::string model_1d, std::string model_3d*/)
{

  Py_NoSiteFlag = 1;

  Py_Initialize();

  PyObject* sysPath = PySys_GetObject((char*)"path");

  PyList_Insert(sysPath,0, PyString_FromString("."));
  PyList_Insert(sysPath,0, PyString_FromString("python27.zip"));
 
  char* argv[]={"minimal.py"};  
  PySys_SetArgv(1,argv); 

  PyObject* PyFileObject = PyFile_FromString(argv[0], "r");
  PyRun_SimpleFile(PyFile_AsFile(PyFileObject), argv[0]);

  Py_Finalize();

  std::cout << "KRATOS TERMINATED CORRECTLY" << std::endl;
  return 0;
}