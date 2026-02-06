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
    py::class_<LinearOperatorType, typename LinearOperatorType::Pointer>(m, "LinearOperator")
        .def(py::init<Parameters>())
        .def("SpMV", &LinearOperatorType::SpMV)
        .def("TransposeSpMV", &LinearOperatorType::TransposeSpMV)
        .def("Clear", &LinearOperatorType::Clear)
        .def("SetNumRows", &LinearOperatorType::SetNumRows)
        .def("SetNumCols", &LinearOperatorType::SetNumCols)
        .def("GetMatrix", py::overload_cast<>(&LinearOperatorType::GetMatrix), py::return_value_policy::reference)
        .def("NumRows", &LinearOperatorType::NumRows)
        .def("NumCols", &LinearOperatorType::NumCols)
        .def("IsMatrixFree", &LinearOperatorType::IsMatrixFree)
        ;

    using SparseMatrixLinearOperatorType = Future::SparseMatrixLinearOperator<Future::SerialLinearAlgebraTraits>;
    py::class_<SparseMatrixLinearOperatorType, typename SparseMatrixLinearOperatorType::Pointer, LinearOperatorType>(m, "SparseMatrixLinearOperator")
        .def(py::init<typename Future::SerialLinearAlgebraTraits::MatrixType&>())
        ;
}

} // namespace Kratos::Future::Python
