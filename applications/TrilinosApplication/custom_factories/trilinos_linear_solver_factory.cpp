//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Vicente Mataix Ferrandiz
//                  Philipp Bucher
//

// Project includes

// Linear solvers
#include "trilinos_linear_solver_factory.h"
#include "linear_solvers/fallback_linear_solver.h"

#ifndef TRILINOS_EXCLUDE_AZTEC_SOLVER
#include "external_includes/aztec_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS_SOLVER
#include "external_includes/amesos_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS2_SOLVER
#include "external_includes/amesos2_solver.h"
#endif

#ifndef TRILINOS_EXCLUDE_ML_SOLVER
#include "external_includes/ml_solver.h"
#endif

#include "external_includes/amgcl_mpi_solver.h"
#include "external_includes/amgcl_mpi_schur_complement_solver.h"
#include "external_includes/trilinos_monotonicity_preserving_solver.h"

namespace Kratos {

void RegisterTrilinosLinearSolvers()
{
    using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

    using TrilinosFallbackLinearSolverType = FallbackLinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
    const static auto TrilinosFallbackLinearSolverFactory = TrilinosLinearSolverFactory<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosFallbackLinearSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("fallback_linear_solver", TrilinosFallbackLinearSolverFactory);

#ifndef TRILINOS_EXCLUDE_AZTEC_SOLVER
    using AztecSolverType = AztecSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType >;
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
    using AmesosSolverType = AmesosSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType >;
    static auto AmesosSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        AmesosSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("amesos",        AmesosSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("klu",           AmesosSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("super_lu_dist", AmesosSolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("mumps",         AmesosSolverFactory);
#endif

#ifndef TRILINOS_EXCLUDE_AMESOS2_SOLVER
    using Amesos2SolverType = Amesos2Solver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType >;
    static auto Amesos2SolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        Amesos2SolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("amesos2",        Amesos2SolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("klu2",           Amesos2SolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("mumps2",         Amesos2SolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("mumsps2",        Amesos2SolverFactory); // NOTE: This is a typo, to be removed, keep it for retrocompatibility
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("super_lu_dist2", Amesos2SolverFactory);
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("basker",         Amesos2SolverFactory);
#endif

#ifndef TRILINOS_EXCLUDE_ML_SOLVER
    using MLSolverType = MultiLevelSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType >;
    static auto MultiLevelSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        MLSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("multi_level", MultiLevelSolverFactory);
#endif

    using AmgclMPISolverType = AmgclMPISolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType >;
    static auto AmgclMPISolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        AmgclMPISolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("amgcl", AmgclMPISolverFactory);

    using AmgclMPISchurComplementSolverType = AmgclMPISchurComplementSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType >;
    static auto AmgclMPISchurComplementSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        AmgclMPISchurComplementSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("amgcl_schur_complement", AmgclMPISchurComplementSolverFactory);

    using TrilinosMonotonicityPreservingSolverType = TrilinosMonotonicityPreservingSolver<TrilinosSparseSpaceType,
        TrilinosLocalSpaceType >;
    static auto TrilinosMonotonicityPreservingSolverFactory = TrilinosLinearSolverFactory<
        TrilinosSparseSpaceType,
        TrilinosLocalSpaceType,
        TrilinosMonotonicityPreservingSolverType>();
    KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER("monotonicity_preserving", TrilinosMonotonicityPreservingSolverFactory);
}

template class KratosComponents<TrilinosLinearSolverFactoryType>;
}
