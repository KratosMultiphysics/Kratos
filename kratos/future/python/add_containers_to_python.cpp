
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "includes/define_python.h"

// Future Extensions
#include "future/python/add_containers_to_python.h"
#include "future/containers/define_linear_algebra_serial.h"
#include "future/containers/linear_system_container.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddContainersToPython(py::module& m)
{
    using LinearSystemContainerType = Future::LinearSystemContainer<Future::SerialLinearAlgebra>;
    py::class_<LinearSystemContainerType, typename LinearSystemContainerType::Pointer>(m, "LinearSystemContainer")
        .def(py::init<>())
        .def("Clear", &LinearSystemContainerType::Clear)
        .def("RequiresEffectiveDofSet", &LinearSystemContainerType::RequiresEffectiveDofSet)
    ;
}

}  // namespace Kratos::Future::Python.

