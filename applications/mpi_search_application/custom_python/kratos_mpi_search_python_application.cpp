//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "mpi_search_application.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos {
namespace Python {

using namespace boost::python;

PYBIND11_MODULE(KratosMPISearchApplication, m)
{

    class_<KratosMPISearchApplication, KratosMPISearchApplication::Pointer,
           KratosApplication>(m, "KratosMPISearchApplication")
           ;

    AddCustomUtilitiesToPython(m);
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(NEIGHBOUR_PARTITION_INDEX)
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
