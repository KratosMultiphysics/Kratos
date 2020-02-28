//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher), based on the work of Jordi Cotela
//

// Project includes
#include "mpi/factories/mpi_linear_solver_factory.h"

// Linear solvers

// Those are in the TrilinosApp! => have to be refactored and moved
// #include "amgcl_mpi_solver.h"
// #include "amgcl_mpi_schur_complement_solver.h"


namespace Kratos {

void RegisterMPILinearSolvers()
{
    // typedef AmgclMPISolver<MPISparseSpaceType,
    //     MPILocalSpaceType > AmgclMPISolverType;
    // static auto AmgclMPISolverFactory = MPILinearSolverFactory<
    //     MPISparseSpaceType,
    //     MPILocalSpaceType,
    //     AmgclMPISolverType>();
    // KRATOS_REGISTER_MPI_LINEAR_SOLVER("amgcl", AmgclMPISolverFactory);

    // typedef AmgclMPISchurComplementSolver<MPISparseSpaceType,
    //     MPILocalSpaceType > AmgclMPISchurComplementSolverType;
    // static auto AmgclMPISchurComplementSolverFactory = MPILinearSolverFactory<
    //     MPISparseSpaceType,
    //     MPILocalSpaceType,
    //     AmgclMPISchurComplementSolverType>();
    // KRATOS_REGISTER_MPI_LINEAR_SOLVER("amgcl_schur_complement", AmgclMPISchurComplementSolverFactory);
}

template class KratosComponents<MPILinearSolverFactoryType>;
}