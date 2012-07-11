#include <Python.h>
#include <iostream>

int main(int argc, char *argv[])
{
  Py_Initialize();
  PySys_SetArgv(argc,argv);

  FILE *fp      = fopen (argv[1],   "r+");

  PyRun_SimpleFile (fp, argv[1]);

  Py_Finalize();
  return 0;
}