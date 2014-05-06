#include <Python.h>
#include <iostream>
#include "stdio.h"
#include "string.h"

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

#if PY_MAJOR_VERSION >= 3

/** Unix-like platform char * to wchar_t conversion. */
#ifdef _WIN32
	#include <windows.h>

    /** Windows char * to wchar_t conversion. */
    wchar_t *nstrws_convert(char *raw)
    {
        int size_needed = MultiByteToWideChar(CP_UTF8, 0, raw, -1, NULL, 0);
        wchar_t *rtn = (wchar_t *) calloc(1, size_needed * sizeof(wchar_t));
        MultiByteToWideChar(CP_UTF8, 0, raw, -1, rtn, size_needed);
        return rtn;
    }
#else
    wchar_t *nstrws_convert(char *raw)
    {
        wchar_t *rtn = (wchar_t *) calloc((strlen(raw) + 1) , sizeof(wchar_t) ); //VALGRIND FINDS A PROBLEM HERE... do not understand why---
        setlocale(LC_ALL,"en_US.UTF-8"); //i don't know if this is really needed... definitely does not look good
        mbstowcs(rtn, raw, strlen(raw));
        return rtn;
    }
#endif

int solve(
	const char * kratos_path,
	const char * model_1d_name, 
	const char * model_3d_name,
	const char * script_name
	)
{
  Py_NoSiteFlag = 1;

  Py_Initialize();

  PyObject* sysPath = PySys_GetObject((char*)"path");

  PyList_Insert(sysPath,0, PyBytes_FromString("."));
  PyList_Insert(sysPath,0, PyBytes_FromString(kratos_path));
  
  char python_lib_path[1024];
  strcpy(python_lib_path, kratos_path);
  strcat(python_lib_path, "/python_embedded_stdlib.zip");
  PyList_Insert(sysPath,0, PyBytes_FromString(python_lib_path));

  char kratos_lib_path[1024];
  strcpy(kratos_lib_path, kratos_path);
  strcat(kratos_lib_path, "/libs");
  PyList_Insert(sysPath,0, PyBytes_FromString(kratos_lib_path));
  
  char* argv[1];
  char script_name_aux[1024];
  strcpy(script_name_aux, script_name);
  argv[0] = script_name_aux;  

  //here we convert argv from char to wchar_t ... this is potentially very problematic
  wchar_t **wchar_argv = (wchar_t **) calloc(1, sizeof(wchar_t *));
  wchar_argv[0] = nstrws_convert(argv[0]);

  PySys_SetArgv(1,wchar_argv);

  int error_code = Py_Main(1,wchar_argv);

  Py_Finalize();

  std::cout << "KRATOS TERMINATED CORRECTLY" << std::endl;
  return 0;
}
#else
int solve(
	const char * kratos_path,
	const char * model_1d_name, 
	const char * model_3d_name,
	const char * script_name
	)
{
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
#endif

#undef KRATOS_DLL_EXPORT

