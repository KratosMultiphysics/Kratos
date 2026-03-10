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
#include "future/preconditioners/preconditioner.h"
#include "future/python/add_preconditioners_to_python.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddPreconditionersToPython(py::module& m)
{
    using PreconditionerType = Future::Preconditioner<Future::SerialLinearAlgebraTraits>;
    py::class_<PreconditionerType, py::smart_holder>(m, "Preconditioner")
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("Initialize", &PreconditionerType::Initialize)
        .def("InitializeSolutionStep", &PreconditionerType::InitializeSolutionStep)
        .def("Apply", &PreconditionerType::Apply)
        .def("ApplyTranspose", &PreconditionerType::ApplyTranspose)
        .def("FinalizeSolutionStep", &PreconditionerType::FinalizeSolutionStep)
        .def("Clear", &PreconditionerType::Clear)
        ;
}

} // namespace Kratos::Future::Python
