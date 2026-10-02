
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

    py::enum_<Future::LinearSystemTags::SparseMatrixTag>(m, "SparseMatrixTag")
        .value("LHS", Future::LinearSystemTags::SparseMatrixTag::LHS)
        .value("MassMatrix", Future::LinearSystemTags::SparseMatrixTag::MassMatrix)
        .value("StiffnessMatrix", Future::LinearSystemTags::SparseMatrixTag::StiffnessMatrix)
        .value("DampingMatrix", Future::LinearSystemTags::SparseMatrixTag::DampingMatrix)
        .value("NumberOfTags", Future::LinearSystemTags::SparseMatrixTag::NumberOfTags)
    ;

    py::enum_<Future::LinearSystemTags::DenseVectorTag>(m, "DenseVectorTag")
        .value("RHS", Future::LinearSystemTags::DenseVectorTag::RHS)
        .value("Dx", Future::LinearSystemTags::DenseVectorTag::Dx)
        .value("Eigvals", Future::LinearSystemTags::DenseVectorTag::Eigvals)
        .value("NumberOfTags", Future::LinearSystemTags::DenseVectorTag::NumberOfTags)
    ;

    py::enum_<Future::LinearSystemTags::DenseMatrixTag>(m, "DenseMatrixTag")
        .value("RHS", Future::LinearSystemTags::DenseMatrixTag::RHS)
        .value("Dx", Future::LinearSystemTags::DenseMatrixTag::Dx)
        .value("Eigvects", Future::LinearSystemTags::DenseMatrixTag::Eigvects)
        .value("NumberOfTags", Future::LinearSystemTags::DenseMatrixTag::NumberOfTags)
    ;

    using LinearSystemType = Future::LinearSystem<Future::SerialLinearAlgebraTraits>;
    py::class_<LinearSystemType, typename LinearSystemType::Pointer>(m, "LinearSystem")
        .def(py::init<>())
        .def(py::init<std::string>(), py::arg("Name"))
        .def(py::init<typename LinearSystemType::MatrixType::Pointer, typename LinearSystemType::VectorType::Pointer, typename LinearSystemType::VectorType::Pointer, std::string>(), py::arg("pLhs"), py::arg("pRhs"), py::arg("pSol"), py::arg("Name"))
        .def(py::init<typename LinearSystemType::LinearOperatorType::UniquePointer, typename LinearSystemType::VectorType::Pointer, typename LinearSystemType::VectorType::Pointer, std::string>(), py::arg("pLinearOperator"), py::arg("pRhs"), py::arg("pSol"), py::arg("Name"))
        .def("SetMatrix", [](LinearSystemType& rLinearSystem, typename Future::SerialLinearAlgebraTraits::MatrixType::Pointer pMatrix, Future::LinearSystemTags::SparseMatrixTag Tag) {
            rLinearSystem.pSetMatrix(pMatrix, Tag); }, py::arg("Matrix"), py::arg("Tag"))
        .def("SetMatrix", [](LinearSystemType& rLinearSystem, typename Kratos::shared_ptr<Future::SerialLinearAlgebraTraits::DenseMatrixType> pMatrix, Future::LinearSystemTags::DenseMatrixTag Tag) {
            rLinearSystem.pSetMatrix(pMatrix, Tag); }, py::arg("Matrix"), py::arg("Tag"))
        .def("SetVector", [](LinearSystemType& rLinearSystem, typename Future::SerialLinearAlgebraTraits::VectorType::Pointer pVector, Future::LinearSystemTags::DenseVectorTag Tag) {
            rLinearSystem.pSetVector(pVector, Tag); }, py::arg("Vector"), py::arg("Tag"))
        .def("SetLinearOperator", [](LinearSystemType& rLinearSystem, typename Future::LinearOperator<Future::SerialLinearAlgebraTraits>::UniquePointer pLinearOperator, Future::LinearSystemTags::SparseMatrixTag Tag) {
            rLinearSystem.pSetLinearOperator(std::move(pLinearOperator), Tag); }, py::arg("LinearOperator"), py::arg("Tag"))
        .def("GetMatrix", [](LinearSystemType& rLinearSystem, Future::LinearSystemTags::SparseMatrixTag Tag) -> typename Future::SerialLinearAlgebraTraits::MatrixType& {
            return *rLinearSystem.pGetMatrix(Tag); }, py::return_value_policy::reference_internal, py::arg("Tag"))
        .def("GetMatrix", [](LinearSystemType& rLinearSystem, Future::LinearSystemTags::DenseMatrixTag Tag) -> typename Future::SerialLinearAlgebraTraits::DenseMatrixType& {
            return *rLinearSystem.pGetMatrix(Tag); }, py::return_value_policy::reference_internal, py::arg("Tag"))
        .def("GetVector", [](LinearSystemType& rLinearSystem, Future::LinearSystemTags::DenseVectorTag Tag) -> typename Future::SerialLinearAlgebraTraits::VectorType& {
            return *rLinearSystem.pGetVector(Tag); }, py::return_value_policy::reference_internal, py::arg("Tag"))
        .def("GetLinearOperator", [](LinearSystemType& rLinearSystem, Future::LinearSystemTags::SparseMatrixTag Tag) -> typename Future::LinearOperator<Future::SerialLinearAlgebraTraits>& {
            return *rLinearSystem.pGetLinearOperator(Tag); }, py::return_value_policy::reference_internal, py::arg("Tag"))
        .def("SetAdditionalData", &LinearSystemType::SetAdditionalData)
        .def_property_readonly("Name", &LinearSystemType::Name)
        .def_property_readonly("HasAdditionalData", &LinearSystemType::HasAdditionalData)
    ;
}

}  // namespace Kratos::Future::Python.

