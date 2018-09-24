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
}

}
}