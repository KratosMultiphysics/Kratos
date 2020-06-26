#if defined(KRATOS_DEBUG) && ! defined(EXCLUDE_EMBEDDED_PYTHON_DEBUG)
  #include <Python.h>
#else
  #ifdef _DEBUG
    #define _DEBUG_DEFINED
    #undef _DEBUG
  #endif
  #include <Python.h>
  #ifdef _DEBUG_DEFINED
    #undef _DEBUG_DEFINED
    #define _DEBUG
  #endif
#endif

#include <iostream>

#if PY_MAJOR_VERSION >= 3
    #include <stdlib.h>
    #include <locale.h>
    #include <string.h>
    #include <wchar.h>


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

    //EXPLORE THE POSSIBILITY OF USING "wmain" instead of "main" in windows.
    //this shall allow avoiding the explicit casting of argv to wchar_argv
    int main(int argc, char *argv[])
    {
        if ( argc < 2)
        {
            std::cout << "MISSING SCRIPT NAME! Usage: " << argv[ 0] << " filename_in.py " <<std::endl;
            return 1;
        }

        //here we convert argv from char to wchar_t ... this is potentially very problematic
        wchar_t **wchar_argv = (wchar_t **) calloc(argc, sizeof(wchar_t *));
        for(int i = 0; i<argc; i++)
        {
                wchar_argv[i] = nstrws_convert(argv[i]);
        }

        Py_NoSiteFlag = 1;
        Py_SetProgramName(wchar_argv[0]);
        Py_Initialize();

        PySys_SetArgv(argc-1, &wchar_argv[1] );

        int error_code = Py_Main(argc,wchar_argv);

        Py_Finalize();

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
        
        //free memory
        //Py_DECREF(sysPath);
        for(int i = 0; i<argc; i++)
            free(wchar_argv[i]);
        free(wchar_argv);
    }
#else
    int main(int argc, char *argv[])
    {
        if ( argc < 2)
        {
            std::cout << "MISSING SCRIPT NAME! Usage: " << argv[ 0] << " filename_in.py " <<std::endl;
            return 1;
        }

        Py_NoSiteFlag = 1;
        Py_SetProgramName(argv[0]);

        Py_Initialize();

        PyObject* sysPath = PySys_GetObject((char*)"path");
        PyList_Insert(sysPath,0, PyString_FromString("."));
        PySys_SetArgv(argc-1,&argv[1]);

        char filename[1024];
        char options[128];
        strcpy(filename,argv[1]);
        strcpy(options,"r");
        PyObject* MyPyFileObject = PyFile_FromString(filename, options);
        int error_code = PyRun_SimpleFile(PyFile_AsFile(MyPyFileObject), argv[1]);

        Py_Finalize();

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
#endif
