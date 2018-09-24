//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Pooyan Dadvand
//

// System includes
#include "includes/define_python.h"
#include "includes/kratos_version.h"

// External includes
#include <pybind11/pybind11.h>

// Module includes
#include "mpi_environment.h"
#include "add_mpi_communicator_to_python.h"

namespace Kratos {
namespace Python {

void module_greet()
{
	std::stringstream header;
	header << "Hello, I am the MPI extension of Kratos Multi-Physics " << KRATOS_VERSION <<" ;-)\n";
    std::cout << header.str();
}

PYBIND11_MODULE(KratosMPI, m)
{
    namespace py = pybind11;

    m.def("Hello",module_greet);

    m.def("MPIInitialize",[](py::list& rPythonArgv) {

        // Converting python sys.argv to c-style argc, argv
        int argc = py::len(rPythonArgv);
        std::cout << argc;
        char** argv = new char*[argc];
        for (int i = 0; i < argc; i++)
        {
            pybind11::str arg_str = py::str(rPythonArgv[i]);
            argv[i] = new char[len(arg_str)];
            std::cout << argv[i];
        }

        MPIEnvironment::Initialize(argc,argv);

        // Deallocating argv
        for (int i = 0; i < argc; i++)
        {
            delete argv[i];
        }
        delete[] argv;
    });

    m.def("MPIFinalize",MPIEnvironment::Finalize);

    AddMPICommunicatorToPython(m);
}

}
}