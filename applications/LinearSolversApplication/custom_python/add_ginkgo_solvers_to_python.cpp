//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "custom_python/add_ginkgo_solvers_to_python.h"
#include "custom_solvers/ginkgo_solver.h"

namespace Kratos::Python {

void AddGinkgoSolversToPython(pybind11::module& m)
{
    using SpaceType = UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SpaceType,  LocalSpaceType>;

    namespace py = pybind11;

    using GinkgoSolverType = GinkgoSolver<SpaceType,  LocalSpaceType>;
    py::class_<GinkgoSolverType, GinkgoSolverType::Pointer, LinearSolverType>(m,"GinkgoSolver")
        .def(py::init<Parameters>())
        .def("__str__", PrintObject<LinearSolverType>)
        .def("CleanSolve", &GinkgoSolverType::CleanSolve)
        ;
}

}  // namespace Kratos::Python.
