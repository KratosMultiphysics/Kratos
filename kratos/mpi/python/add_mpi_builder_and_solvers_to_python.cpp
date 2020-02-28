
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
#include "mpi/python/add_mpi_builder_and_solvers_to_python.h"

#include "mpi/spaces/amgcl_mpi_space.h"
#include "spaces/ublas_space.h"

#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "linear_solvers/linear_solver.h"

// BuilderAndSolvers
#include "mpi/solving_strategies/builder_and_solvers/mpi_block_builder_and_solver.h"


namespace Kratos {
namespace Python {

void AddMPIBuilderAndSolversToPython(pybind11::module& m)
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
    typedef BuilderAndSolver< MPISparseSpaceType, MPILocalSpaceType, MPILinearSolverType > MPIBuilderAndSolverType;

    // Base BuilderAndSolver
    py::class_< MPIBuilderAndSolverType, typename MPIBuilderAndSolverType::Pointer >(m, "MPIResidualBasedBuilderAndSolver")
        .def(py::init<MPILinearSolverType::Pointer> () )
        .def( "SetCalculateReactionsFlag", &MPIBuilderAndSolverType::SetCalculateReactionsFlag )
        .def( "GetCalculateReactionsFlag", &MPIBuilderAndSolverType::GetCalculateReactionsFlag )
        .def( "SetDofSetIsInitializedFlag", &MPIBuilderAndSolverType::SetDofSetIsInitializedFlag )
        .def( "GetDofSetIsInitializedFlag", &MPIBuilderAndSolverType::GetDofSetIsInitializedFlag )
        .def( "SetReshapeMatrixFlag", &MPIBuilderAndSolverType::SetReshapeMatrixFlag )
        .def( "GetReshapeMatrixFlag", &MPIBuilderAndSolverType::GetReshapeMatrixFlag )
        .def( "GetEquationSystemSize", &MPIBuilderAndSolverType::GetEquationSystemSize )
        .def( "BuildLHS", &MPIBuilderAndSolverType::BuildLHS )
        .def( "BuildRHS", &MPIBuilderAndSolverType::BuildRHS )
        .def( "Build", &MPIBuilderAndSolverType::Build )
        .def( "SystemSolve", &MPIBuilderAndSolverType::SystemSolve )
        .def( "BuildAndSolve", &MPIBuilderAndSolverType::BuildAndSolve )
        .def( "BuildRHSAndSolve", &MPIBuilderAndSolverType::BuildRHSAndSolve )
        .def( "ApplyDirichletConditions", &MPIBuilderAndSolverType::ApplyDirichletConditions )
        .def( "SetUpDofSet", &MPIBuilderAndSolverType::SetUpDofSet )
        .def( "GetDofSet", &MPIBuilderAndSolverType::GetDofSet, py::return_value_policy::reference_internal )
        .def( "SetUpSystem", &MPIBuilderAndSolverType::SetUpSystem )
        .def( "ResizeAndInitializeVectors", &MPIBuilderAndSolverType::ResizeAndInitializeVectors )
        .def( "InitializeSolutionStep", &MPIBuilderAndSolverType::InitializeSolutionStep )
        .def( "FinalizeSolutionStep", &MPIBuilderAndSolverType::FinalizeSolutionStep )
        .def( "CalculateReactions", &MPIBuilderAndSolverType::CalculateReactions )
        .def( "Clear", &MPIBuilderAndSolverType::Clear )
        .def( "SetEchoLevel", &MPIBuilderAndSolverType::SetEchoLevel )
        .def( "GetEchoLevel", &MPIBuilderAndSolverType::GetEchoLevel )
        ;

    // BuilderAndSolvers
    typedef MPIBlockBuilderAndSolver< MPISparseSpaceType, MPILocalSpaceType, MPILinearSolverType > MPIBlockBuilderAndSolverType;
    py::class_< MPIBlockBuilderAndSolverType, typename MPIBlockBuilderAndSolverType::Pointer,MPIBuilderAndSolverType >(m, "MPIBlockBuilderAndSolver")
        .def(py::init<MPILinearSolverType::Pointer> () )
        .def(py::init<MPILinearSolverType::Pointer, Parameters> () )
        ;
}

} // namespace Python.
} // Namespace Kratos
