
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:     Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "mpi/python/add_mpi_solving_strategies_to_python.h"

#include "mpi/spaces/amgcl_mpi_space.h"
#include "spaces/ublas_space.h"

#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "linear_solvers/linear_solver.h"

// SolvingStrategies
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"


namespace Kratos {
namespace Python {

void AddMPISolvingStrategiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Amgcl type definitions
    typedef amgcl::backend::builtin<double> Backend;
    typedef amgcl::mpi::distributed_matrix<Backend> amgcl_mpi_matrix;
    typedef typename Backend::vector amgcl_mpi_vector;

    // Type definitions
    typedef AmgclMPISpace<amgcl_mpi_matrix, amgcl_mpi_vector> MPISparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> MPILocalSpaceType;
    typedef LinearSolver<MPISparseSpaceType, MPILocalSpaceType > MPILinearSolverType;
    typedef ConvergenceCriteria< MPISparseSpaceType, MPILocalSpaceType > MPIConvergenceCriteria;

    typedef SolvingStrategy< MPISparseSpaceType, MPILocalSpaceType, MPILinearSolverType > MPIBaseSolvingStrategyType;
    typedef Scheme< MPISparseSpaceType, MPILocalSpaceType > MPISchemeType;
    typedef BuilderAndSolver< MPISparseSpaceType, MPILocalSpaceType, MPILinearSolverType > MPIBuilderAndSolverType;

    // Base SolvingStrategy
    py::class_< MPIBaseSolvingStrategyType, typename MPIBaseSolvingStrategyType::Pointer >(m, "MPISolvingStrategy")
        .def(py::init< ModelPart&, bool >())
        .def("Initialize", &MPIBaseSolvingStrategyType::Initialize)
        .def("Check", &MPIBaseSolvingStrategyType::Check)
        .def("Clear", &MPIBaseSolvingStrategyType::Clear)
        .def("InitializeSolutionStep", &MPIBaseSolvingStrategyType::InitializeSolutionStep)
        .def("Predict", &MPIBaseSolvingStrategyType::Predict)
        .def("SolveSolutionStep", &MPIBaseSolvingStrategyType::SolveSolutionStep)
        .def("FinalizeSolutionStep", &MPIBaseSolvingStrategyType::FinalizeSolutionStep)
        .def("IsConverged", &MPIBaseSolvingStrategyType::IsConverged)
        .def("SetEchoLevel", &MPIBaseSolvingStrategyType::SetEchoLevel)
        .def("GetEchoLevel", &MPIBaseSolvingStrategyType::GetEchoLevel)
        .def("SetRebuildLevel", &MPIBaseSolvingStrategyType::SetRebuildLevel)
        .def("GetRebuildLevel", &MPIBaseSolvingStrategyType::GetRebuildLevel)
        .def("SetMoveMeshFlag", &MPIBaseSolvingStrategyType::SetMoveMeshFlag)
        .def("MoveMeshFlag", &MPIBaseSolvingStrategyType::MoveMeshFlag)
        .def("MoveMesh", &MPIBaseSolvingStrategyType::MoveMesh)
        ;

    // SolvingStrategies
    typedef ResidualBasedLinearStrategy< MPISparseSpaceType, MPILocalSpaceType, MPILinearSolverType> MPILinearStrategy;
    py::class_< MPILinearStrategy , typename MPILinearStrategy::Pointer, MPIBaseSolvingStrategyType >
    (m,"MPILinearStrategy")
        .def(py::init< ModelPart&, MPISchemeType::Pointer, MPILinearSolverType::Pointer, MPIBuilderAndSolverType::Pointer, bool, bool, bool, bool >());

    typedef ResidualBasedNewtonRaphsonStrategy< MPISparseSpaceType, MPILocalSpaceType, MPILinearSolverType> MPINewtonRaphsonStrategy;
    py::class_< MPINewtonRaphsonStrategy , typename MPINewtonRaphsonStrategy::Pointer, MPIBaseSolvingStrategyType >
    (m,"MPINewtonRaphsonStrategy")
        .def(py::init< ModelPart&, MPISchemeType::Pointer, MPILinearSolverType::Pointer, MPIConvergenceCriteria::Pointer, MPIBuilderAndSolverType::Pointer, int, bool, bool, bool >());
}

} // namespace Python.
} // Namespace Kratos
