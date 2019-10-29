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

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rans_variable_utilities.h"

namespace Kratos
{
namespace Python
{
void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def_submodule("RansVariableUtilities")
        .def("ClipScalarVariable", &RansVariableUtilities::ClipScalarVariable)
        .def("GetMinimumScalarValue", &RansVariableUtilities::GetMinimumScalarValue)
        .def("GetMaximumScalarValue", &RansVariableUtilities::GetMaximumScalarValue)
        .def("CopyNodalSolutionStepVariablesList",
             &RansVariableUtilities::CopyNodalSolutionStepVariablesList);
}

} // namespace Python.
} // Namespace Kratos
