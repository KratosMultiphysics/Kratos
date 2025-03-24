
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
#include "includes/define_python.h"

// Future Extensions
#include "future/python/add_linear_solvers_to_python.h"
#include "future/linear_solvers/linear_solver.h"
#include "future/linear_solvers/iterative_solver.h"
#include "future/linear_solvers/skyline_lu_factorization_solver.h"

namespace Kratos::Future::Python
{

namespace py = pybind11;

void AddLinearSolversToPython(py::module& m)
{

    using LinearSolverType = Future::LinearSolver<>;
    py::class_<LinearSolverType, typename LinearSolverType::Pointer>(m, "LinearSolver")
        .def(py::init<>())
        .def("Initialize", &LinearSolverType::Initialize)
        // .def("__str__", PrintObject<LinearSolverType>)
        .def("GetIterationsNumber", &LinearSolverType::GetIterationsNumber)
    ;

    using DirectSolverType = Future::DirectSolver<>;
    py::class_<DirectSolverType, typename DirectSolverType::Pointer, LinearSolverType>(m, "DirectSolver")
        .def(py::init<>())
        .def(py::init<Parameters>())
        // .def("__str__", PrintObject<DirectSolverType>)
    ;

    using SkylineLUFactorizationSolverType = Future::SkylineLUFactorizationSolver<>;
    py::class_<SkylineLUFactorizationSolverType, SkylineLUFactorizationSolverType::Pointer, DirectSolverType>(m, "SkylineLUFactorizationSolver")
        .def(py::init<>())
        .def(py::init<Parameters>())
        // .def("__str__", PrintObject<SkylineLUFactorizationSolverType>)
    ;
}

}  // namespace Kratos::Future::Python.

