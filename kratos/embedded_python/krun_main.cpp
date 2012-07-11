#include <Python.h>
#include <iostream>

int main(int argc, char *argv[])
{
	if ( argc < 2) {
		std::cout << "MISSING SCRIPT NAME! Usage: " << argv[ 0] << " filename_in.py " <<std::endl;
		return 1;
	}
  Py_Initialize();
  
  PySys_SetArgv(argc,argv);

 // FILE *fp      = fopen (argv[1],   "r");
 // if ( !fp) {
	 // std::cout << "Error opening file: " << argv[ 1] << "\n";
	 // return 1;
  //}

  PyObject* PyFileObject = PyFile_FromString(argv[1], "r");
  PyRun_SimpleFile(PyFile_AsFile(PyFileObject), argv[1]);
  //PyRun_SimpleFile (fp, argv[1]);

  Py_Finalize();
  //fclose( fp);

  std::cout << "KRATOS TERMINATED CORRECTLY" << std::endl;
  return 0;
}