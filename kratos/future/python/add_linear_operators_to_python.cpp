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

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/kratos_parameters.h"

// Future Extensions
#include "future/containers/define_linear_algebra_serial.h"
#include "future/linear_operators/linear_operator.h"
#include "future/linear_operators/sparse_matrix_linear_operator.h"
#include "future/python/add_linear_operators_to_python.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddLinearOperatorsToPython(py::module& m)
{
    using LinearOperatorType = Future::LinearOperator<Future::SerialLinearAlgebraTraits>;
    py::class_<LinearOperatorType, py::smart_holder>(m, "LinearOperator")
        .def(py::init<Parameters>())
        .def(py::init<const std::pair<std::size_t, std::size_t>>())
        .def("SpMV", &LinearOperatorType::SpMV)
        .def("TransposeSpMV", &LinearOperatorType::TransposeSpMV)
        .def("GetMatrix", py::overload_cast<>(&LinearOperatorType::GetMatrix), py::return_value_policy::reference_internal)
        .def_property_readonly("Shape", &LinearOperatorType::Shape)
        .def_property_readonly("Size1", &LinearOperatorType::Size1)
        .def_property_readonly("Size2", &LinearOperatorType::Size2)
        .def_property_readonly("IsMatrixFree", &LinearOperatorType::IsMatrixFree)
        ;

    using SparseMatrixLinearOperatorType = Future::SparseMatrixLinearOperator<Future::SerialLinearAlgebraTraits>;
    py::class_<SparseMatrixLinearOperatorType, py::smart_holder, LinearOperatorType>(m, "SparseMatrixLinearOperator")
        .def(py::init<typename Future::SerialLinearAlgebraTraits::MatrixType&>())
        ;
}

} // namespace Kratos::Future::Python
