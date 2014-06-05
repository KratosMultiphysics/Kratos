#include <Python.h>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include "string.h"
#include <cstdio>

//#include "cvp.h"
#ifdef _WINDLL
#define KRATOS_DLL_EXPORT __declspec(dllexport)
#else
#define KRATOS_DLL_EXPORT
#endif
KRATOS_DLL_EXPORT int solve(
	const char * kratos_path,
	const char * model_1d_name, 
	const char * model_3d_name,
	const char * script_name
	);

int solve(
	const char * kratos_path,
	const char * model_1d_name, 
	const char * model_3d_name,
	const char * script_name
	)
{
//std::ofstream file("logfile.out");
//std::cout.rdbuf(file.rdbuf());

dup2(1, 2);

freopen("output.txt","w",stdout);


//freopen("outerr.txt","w",stderr);


std::cout << "entered in function" << std::endl;
//freopen("output_err.txt","w",stderr);

  Py_NoSiteFlag = 1;

  Py_Initialize();

  PyObject* sysPath = PySys_GetObject((char*)"path");

  PyList_Insert(sysPath,0, PyString_FromString("."));
  PyList_Insert(sysPath,0, PyString_FromString(kratos_path));
  
  char python_lib_path[1024];
  strcpy(python_lib_path, kratos_path);
  strcat(python_lib_path, "/python27.zip");
  PyList_Insert(sysPath,0, PyString_FromString(python_lib_path));

  char kratos_lib_path[1024];
  strcpy(kratos_lib_path, kratos_path);
  strcat(kratos_lib_path, "/libs");
  PyList_Insert(sysPath,0, PyString_FromString(kratos_lib_path));
  
  char* argv[1];
  char script_name_aux[1024];
  strcpy(script_name_aux, script_name);
  argv[0] = script_name_aux;  
  PySys_SetArgv(1,argv);


  PyObject* PyFileObject = PyFile_FromString(argv[0], "r");
  PyRun_SimpleFile(PyFile_AsFile(PyFileObject), argv[0]);

  Py_Finalize();

  std::cout << "KRATOS TERMINATED CORRECTLY" << std::endl;
  return 0;
}

#undef KRATOS_DLL_EXPORT

