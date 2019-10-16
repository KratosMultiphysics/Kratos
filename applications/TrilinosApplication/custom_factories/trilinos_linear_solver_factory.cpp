//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:   Vicente Mataix Ferrandiz
//                  Philipp Bucher
//


// // Project includes

// // Linear solvers
#include "trilinos_linear_solver_factory.h"

#ifndef TRILINOS_EXCLUDE_AZTEC_SOLVER
#include "external_includes/aztec_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS_SOLVER
#include "external_includes/amesos_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_ML_SOLVER
#include "external_includes/ml_solver.h"
#endif

#include "external_includes/amgcl_mpi_solver.h"
#include "external_includes/amgcl_mpi_schur_complement_solver.h"

namespace Kratos {

void RegisterTrilinosLinearSolvers()
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

#ifndef TRILINOS_EXCLUDE_AZTEC_SOLVER
    typedef AztecSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType > AztecSolverType;
    static auto AztecSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        AztecSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("aztec",    AztecSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("cg",       AztecSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("bicgstab", AztecSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("gmres",    AztecSolverFactory);
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS_SOLVER
    typedef AmesosSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType > AmesosSolverType;
    static auto AmesosSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        AmesosSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("amesos",        AmesosSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("klu",           AmesosSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("super_lu_dist", AmesosSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("mumps",         AmesosSolverFactory);
#endif

#ifndef TRILINOS_EXCLUDE_ML_SOLVER
    typedef MultiLevelSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType > MLSolverType;
    static auto MultiLevelSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        MLSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("multi_level", MultiLevelSolverFactory);
#endif

    typedef AmgclMPISolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType > AmgclMPISolverType;
    static auto AmgclMPISolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        AmgclMPISolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("amgcl", AmgclMPISolverFactory);

    typedef AmgclMPISchurComplementSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType > AmgclMPISchurComplementSolverType;
    static auto AmgclMPISchurComplementSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        AmgclMPISchurComplementSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("amgcl_schur_complement", AmgclMPISchurComplementSolverFactory);
}

template class KratosComponents<TrilinosLinearSolverFactoryType>;
}