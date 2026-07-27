//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "meshioplusplus_application.h"
#include "custom_python/add_custom_io_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos::Python
{

PYBIND11_MODULE(KratosMeshioPlusPlusApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosMeshioPlusPlusApplication, KratosMeshioPlusPlusApplication::Pointer, KratosApplication>(
        m, "KratosMeshioPlusPlusApplication")
        .def(py::init<>())
        ;

    AddCustomIOToPython(m);
    AddCustomUtilitiesToPython(m);
}

} // namespace Kratos::Python

#endif // KRATOS_PYTHON defined
