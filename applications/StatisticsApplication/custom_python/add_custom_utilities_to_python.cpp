//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
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
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<int>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<double>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<array_1d<double, 3>>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<Vector>)
        .def("GetNormMethod", &MethodUtilities::GetNormMethod<Matrix>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<int>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<double>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<array_1d<double, 3>>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<Vector>)
        .def("RaiseToPower", &MethodUtilities::RaiseToPower<Matrix>)
        ;


}

} // namespace Python.
} // Namespace Kratos
