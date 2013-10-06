#include <Python.h>
#include <iostream>

int main(int argc, char *argv[])
{
	if ( argc < 2) {
		std::cout << "MISSING SCRIPT NAME! Usage: " << argv[ 0] << " filename_in.py " <<std::endl;
		return 1;
	}

  Py_NoSiteFlag = 1;

  Py_Initialize();

  PyObject* sysPath = PySys_GetObject((char*)"path");
//   PyList_Append(sysPath, PyString_FromString("."));
//  // PyList_Append(sysPath, PyString_FromString("libs"));
//   PyList_Append(sysPath, PyString_FromString("python27.zip"));
//   //PySys_SetPath("python27.zip:.");
//   //PySys_SetPath(".");
PyList_Insert(sysPath,0, PyString_FromString("."));
PyList_Insert(sysPath,0, PyString_FromString("python27.zip"));
  PySys_SetArgv(argc,argv);

 // FILE *fp      = fopen (argv[1],   "r");
 // if ( !fp) {
	 // std::cout << "Error opening file: " << argv[ 1] << "\n";
	 // return 1;
  //}

  PyObject* PyFileObject = PyFile_FromString(argv[1], "r");
  int error_code = PyRun_SimpleFile(PyFile_AsFile(PyFileObject), argv[1]);
  
  //PyRun_SimpleFile (fp, argv[1]);

  Py_Finalize();
  //fclose( fp);

  if(error_code != 0) 
  {
	  std::cout << "KRATOS TERMINATED WITH ERROR" << std::endl;
	  return 1;
  }
  else
  {
		std::cout << "KRATOS TERMINATED CORRECTLY" << std::endl;
		return 0;
  }
}