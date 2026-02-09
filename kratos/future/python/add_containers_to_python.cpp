
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
#include "future/containers/eigenvalue_system.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddContainersToPython(py::module& m)
{

    py::enum_<Future::LinearSystemTag>(m, "LinearSystemTag")
        .value("LHS", Future::LinearSystemTag::LHS)
        .value("RHS", Future::LinearSystemTag::RHS)
        .value("Dx", Future::LinearSystemTag::Dx)
        .value("M", Future::LinearSystemTag::M)
        .value("K", Future::LinearSystemTag::K)
        .value("C", Future::LinearSystemTag::C)
        .value("Eigvals", Future::LinearSystemTag::Eigvals)
        .value("Eigvects", Future::LinearSystemTag::Eigvects)
        .value("NumberOfTags", Future::LinearSystemTag::NumberOfTags)
    ;

    using LinearSystemType = Future::LinearSystem<Future::SerialLinearAlgebraTraits>;
    py::class_<LinearSystemType, typename LinearSystemType::Pointer>(m, "LinearSystem")
        .def(py::init<>())
        .def(py::init<std::string>(), py::arg("Name"))
        .def(py::init<typename LinearSystemType::MatrixType::Pointer, typename LinearSystemType::VectorType::Pointer, typename LinearSystemType::VectorType::Pointer, std::string>(), py::arg("pLhs"), py::arg("pRhs"), py::arg("pSol"), py::arg("Name"))
        .def(py::init<typename LinearSystemType::LinearOperatorType::UniquePointer, typename LinearSystemType::VectorType::Pointer, typename LinearSystemType::VectorType::Pointer, std::string>(), py::arg("pLinearOperator"), py::arg("pRhs"), py::arg("pSol"), py::arg("Name"))
        .def("GetMatrix", [](LinearSystemType& rLinearSystem, Future::LinearSystemTag Tag) -> typename Future::SerialLinearAlgebraTraits::MatrixType& { return rLinearSystem.GetMatrix(Tag); }, py::return_value_policy::reference_internal)
        .def("GetVector", [](LinearSystemType& rLinearSystem, Future::LinearSystemTag Tag) -> typename Future::SerialLinearAlgebraTraits::VectorType& { return rLinearSystem.GetVector(Tag); }, py::return_value_policy::reference_internal)
        .def("GetLinearOperator", [](LinearSystemType& rLinearSystem, Future::LinearSystemTag Tag) -> typename Future::LinearOperator<Future::SerialLinearAlgebraTraits>& { return rLinearSystem.GetLinearOperator(Tag); }, py::return_value_policy::reference_internal)
        .def("SetAdditionalData", &LinearSystemType::SetAdditionalData)
        .def_property_readonly("Name", &LinearSystemType::Name)
        .def_property_readonly("HasAdditionalData", &LinearSystemType::HasAdditionalData)
    ;

    using EigenvalueSystemType = Future::EigenvalueSystem<Future::SerialLinearAlgebraTraits>;
    py::class_<EigenvalueSystemType, typename EigenvalueSystemType::Pointer>(m, "EigenvalueSystem")
        .def(py::init<>())
        .def(py::init<typename EigenvalueSystemType::MatrixType::Pointer, typename EigenvalueSystemType::MatrixType::Pointer, typename EigenvalueSystemType::VectorType::Pointer, typename EigenvalueSystemType::DenseMatrixPointerType, std::string>(), py::arg("pK"), py::arg("pM"), py::arg("pEigenvalues"), py::arg("pEigenvectors"), py::arg("Name") = "")
        .def("GetStiffnessMatrix", [](EigenvalueSystemType& rEigenvalueSystem)->typename Future::SerialLinearAlgebraTraits::MatrixType& { return rEigenvalueSystem.GetStiffnessMatrix(); }, py::return_value_policy::reference)
        .def("GetMassMatrix", [](EigenvalueSystemType& rEigenvalueSystem)->typename Future::SerialLinearAlgebraTraits::MatrixType& { return rEigenvalueSystem.GetMassMatrix(); }, py::return_value_policy::reference)
        .def("GetEigenvalues", [](EigenvalueSystemType& rEigenvalueSystem)->typename Future::SerialLinearAlgebraTraits::VectorType& { return rEigenvalueSystem.GetEigenvalues(); }, py::return_value_policy::reference)
        .def("GetEigenvectors", [](EigenvalueSystemType& rEigenvalueSystem)->typename Future::SerialLinearAlgebraTraits::DenseMatrixType& { return rEigenvalueSystem.GetEigenvectors(); }, py::return_value_policy::reference)
        .def("Name", &EigenvalueSystemType::Name)
    ;

}

}  // namespace Kratos::Future::Python.

