//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "mesh_moving_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosMeshMovingApplication,m) {
    namespace py = pybind11;

    py::class_<KratosMeshMovingApplication, KratosMeshMovingApplication::Pointer,
        KratosApplication>(m,"KratosMeshMovingApplication")
        .def(py::init<>());

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
