//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "add_amgcl_solver_to_python.h"
#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/amgcl_ns_solver.h"


namespace Kratos::Python {


void  AddAMGCLSolverToPython(pybind11::module& m)
{
    using SpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SpaceType, LocalSpaceType>;

    namespace py = pybind11;

    using AMGCLSolverType = AMGCLSolver<SpaceType, LocalSpaceType>;
    py::class_<AMGCLSolverType, std::shared_ptr<AMGCLSolverType>, LinearSolverType>(m, "AMGCLSolver")
        .def(py::init<const std::string&, const std::string&, double, int, int, int>(),
             py::arg("smoother_name"),
             py::arg("solver_name"),
             py::arg("tolerance"),
             py::arg("max_iterations"),
             py::arg("verbosity"),
             py::arg("gmres_size") = 50)
        .def(py::init<const std::string&, const std::string&, const std::string&, double, int, int, int, bool>(),
             py::arg("smoother_name"),
             py::arg("solver_name"),
             py::arg("coarsening_name"),
             py::arg("tolerance"),
             py::arg("max_iterations"),
             py::arg("verbosity"),
             py::arg("gmres_size") = 50,
             py::arg("provide_coordinates") = false)
        .def(py::init<>())
        .def(py::init<Parameters>(),
             py::arg("settings"))
        ;

    using AMGCL_NS_SolverType = AMGCL_NS_Solver<SpaceType, LocalSpaceType>;
    py::class_<AMGCL_NS_SolverType, std::shared_ptr<AMGCL_NS_SolverType>, LinearSolverType>(m, "AMGCL_NS_Solver")
        .def(py::init<Parameters>())
        ;
}

}  // namespace Kratos::Python.
