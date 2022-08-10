//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_processes/max_singular_value_decomposition_process.h"

// Include base h
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MaxSingularValueDecompositionProcess, MaxSingularValueDecompositionProcess::Pointer, Process>(m, "MaxSingularValueDecompositionProcess")
        .def(py::init<Model&, Parameters&>());

}
} // namespace Python
} // namespace Kratos
