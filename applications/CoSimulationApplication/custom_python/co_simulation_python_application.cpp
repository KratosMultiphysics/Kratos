//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "co_simulation_application.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosCoSimulationApplication,m)
{
    pybind11::class_<KratosCoSimulationApplication,
        KratosCoSimulationApplication::Pointer,
        KratosApplication >(m,"KratosCoSimulationApplication")
        .def(pybind11::init<>())
        ;

}
    // We can add a second module for stand alone users.
}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
