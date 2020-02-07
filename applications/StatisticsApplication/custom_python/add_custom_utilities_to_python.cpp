//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

// Application includes
#include "custom_utilities/method_utilities.h"

namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def_submodule("MethodUtilities")
        .def("GetNormMethod", &MethodsUtilities::GetNormMethod<int>)
        .def("GetNormMethod", &MethodsUtilities::GetNormMethod<double>)
        .def("GetNormMethod", &MethodsUtilities::GetNormMethod<array_1d<double, 3>>)
        .def("GetNormMethod", &MethodsUtilities::GetNormMethod<Vector>)
        .def("GetNormMethod", &MethodsUtilities::GetNormMethod<Matrix>)
        ;


}

} // namespace Python.
} // Namespace Kratos
