
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
#include "future/containers/implicit_strategy_data_container.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddContainersToPython(py::module& m)
{
    using ImplicitStrategyDataContainerType = Future::ImplicitStrategyDataContainer<Future::SerialLinearAlgebraTraits>;
    py::class_<ImplicitStrategyDataContainerType, typename ImplicitStrategyDataContainerType::Pointer>(m, "ImplicitStrategyDataContainer")
        .def(py::init<>())
        .def("Clear", &ImplicitStrategyDataContainerType::Clear)
        .def("RequiresEffectiveDofSet", &ImplicitStrategyDataContainerType::RequiresEffectiveDofSet)
    ;
}

}  // namespace Kratos::Future::Python.

