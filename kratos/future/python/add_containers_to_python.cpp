
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
#include "future/containers/linear_system.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddContainersToPython(py::module& m)
{

    py::enum_<Future::SparseMatrixTag>(m, "SparseMatrixTag")
        .value("LHS", Future::SparseMatrixTag::LHS)
        .value("MassMatrix", Future::SparseMatrixTag::MassMatrix)
        .value("StiffnessMatrix", Future::SparseMatrixTag::StiffnessMatrix)
        .value("DampingMatrix", Future::SparseMatrixTag::DampingMatrix)
        .value("NumberOfTags", Future::SparseMatrixTag::NumberOfTags)
    ;

    py::enum_<Future::DenseVectorTag>(m, "DenseVectorTag")
        .value("RHS", Future::DenseVectorTag::RHS)
        .value("Dx", Future::DenseVectorTag::Dx)
        .value("Eigvals", Future::DenseVectorTag::Eigvals)
        .value("NumberOfTags", Future::DenseVectorTag::NumberOfTags)
    ;

    py::enum_<Future::DenseMatrixTag>(m, "DenseMatrixTag")
        .value("RHS", Future::DenseMatrixTag::RHS)
        .value("Dx", Future::DenseMatrixTag::Dx)
        .value("Eigvects", Future::DenseMatrixTag::Eigvects)
        .value("NumberOfTags", Future::DenseMatrixTag::NumberOfTags)
    ;

    using LinearSystemType = Future::LinearSystem<Future::SerialLinearAlgebraTraits>;
    py::class_<LinearSystemType, typename LinearSystemType::Pointer>(m, "LinearSystem")
        .def(py::init<>())
        .def(py::init<std::string>(), py::arg("Name"))
        .def(py::init<typename LinearSystemType::MatrixType::Pointer, typename LinearSystemType::VectorType::Pointer, typename LinearSystemType::VectorType::Pointer, std::string>(), py::arg("pLhs"), py::arg("pRhs"), py::arg("pSol"), py::arg("Name"))
        .def(py::init<typename LinearSystemType::LinearOperatorType::UniquePointer, typename LinearSystemType::VectorType::Pointer, typename LinearSystemType::VectorType::Pointer, std::string>(), py::arg("pLinearOperator"), py::arg("pRhs"), py::arg("pSol"), py::arg("Name"))
        .def("GetMatrix", [](LinearSystemType& rLinearSystem, Future::SparseMatrixTag Tag) -> typename Future::SerialLinearAlgebraTraits::MatrixType& { return rLinearSystem.GetMatrix(Tag); }, py::return_value_policy::reference_internal)
        .def("GetMatrix", [](LinearSystemType& rLinearSystem, Future::DenseMatrixTag Tag) -> typename Future::SerialLinearAlgebraTraits::DenseMatrixType& { return rLinearSystem.GetMatrix(Tag); }, py::return_value_policy::reference_internal)
        .def("GetVector", [](LinearSystemType& rLinearSystem, Future::DenseVectorTag Tag) -> typename Future::SerialLinearAlgebraTraits::VectorType& { return rLinearSystem.GetVector(Tag); }, py::return_value_policy::reference_internal)
        .def("GetLinearOperator", [](LinearSystemType& rLinearSystem, Future::SparseMatrixTag Tag) -> typename Future::LinearOperator<Future::SerialLinearAlgebraTraits>& { return rLinearSystem.GetLinearOperator(Tag); }, py::return_value_policy::reference_internal)
        .def("SetAdditionalData", &LinearSystemType::SetAdditionalData)
        .def_property_readonly("Name", &LinearSystemType::Name)
        .def_property_readonly("HasAdditionalData", &LinearSystemType::HasAdditionalData)
    ;
}

}  // namespace Kratos::Future::Python.

