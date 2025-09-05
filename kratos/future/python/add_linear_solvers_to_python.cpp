
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
#include "future/linear_solvers/amgcl_solver.h"
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

    using SkylineLUFactorizationSolverType = Future::SkylineLUFactorizationSolver<CsrMatrix<>, SystemVector<>>;
    py::class_<SkylineLUFactorizationSolverType, SkylineLUFactorizationSolverType::Pointer, DirectSolverType>(m, "SkylineLUFactorizationSolver")
        .def(py::init<>())
        .def(py::init<Parameters>())
        // .def("__str__", PrintObject<SkylineLUFactorizationSolverType>)
    ;

    py::enum_<Future::AMGCLSmoother>(m, "AMGCLSmoother")
        .value("SPAI0", SPAI0)
        .value("SPAI1", SPAI1)
        .value("ILU0", ILU0)
        .value("DAMPED_JACOBI", DAMPED_JACOBI)
        .value("GAUSS_SEIDEL", GAUSS_SEIDEL)
        .value("CHEBYSHEV", CHEBYSHEV)
    ;

    py::enum_<Future::AMGCLIterativeSolverType>(m, "AMGCLIterativeSolverType")
        .value("GMRES", GMRES)
        .value("LGMRES", LGMRES)
        .value("FGMRES", FGMRES)
        .value("BICGSTAB", BICGSTAB)
        .value("CG", CG)
        .value("BICGSTAB_WITH_GMRES_FALLBACK", BICGSTAB_WITH_GMRES_FALLBACK)
        .value("BICGSTAB2", BICGSTAB2)
    ;

    py::enum_<Future::AMGCLCoarseningType>(m, "AMGCLCoarseningType")
        .value("RUGE_STUBEN", RUGE_STUBEN)
        .value("AGGREGATION", AGGREGATION)
        .value("SA", SA)
        .value("SA_EMIN", SA_EMIN)
    ;

    using AMGCLSolverType = Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>;
    py::class_<AMGCLSolverType, typename AMGCLSolverType::Pointer, LinearSolverType>(m, "AMGCLSolver")
        .def(py::init<AMGCLSmoother, AMGCLIterativeSolverType, double, int, int, int>())
        .def(py::init<AMGCLSmoother, AMGCLIterativeSolverType, AMGCLCoarseningType, double, int, int, int, bool>())
        .def(py::init<>())
        .def(py::init<Parameters>())
        .def("GetResidualNorm", &AMGCLSolverType::GetResidualNorm)
    ;

}

}  // namespace Kratos::Future::Python.

