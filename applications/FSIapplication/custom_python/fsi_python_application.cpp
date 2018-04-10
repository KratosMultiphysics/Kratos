//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "fsi_application.h"
#include "custom_python/add_convergence_accelerators_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_mappers_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

PYBIND11_MODULE(KratosFSIApplication,m)
{

    class_<KratosFSIApplication,
           KratosFSIApplication::Pointer,
           KratosApplication>(m,"KratosFSIApplication")
           ;

    AddCustomUtilitiesToPython(m);
    AddMappersToPython(m);
    AddConvergenceAcceleratorsToPython(m);

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
