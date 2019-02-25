//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "mixedelement_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosMixedElementApplication,m)
{

    py::class_<KratosMixedElementApplication,
           KratosMixedElementApplication::Pointer,
           KratosApplication > (m,"KratosMixedElementApplication")
           ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SZ)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SXY)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SXZ)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,SYZ)

}


} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
