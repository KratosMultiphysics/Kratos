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
#include "includes/define.h"
#include "factories/standard_linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/amgcl_ns_solver.h"
#include "linear_solvers/scaling_solver.h"
#include "linear_solvers/monotonicity_preserving_solver.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    void RegisterLinearSolvers()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;
        typedef TUblasSparseSpace<std::complex<double>> ComplexSpaceType;
        typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSpaceType;
//         typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
//         typedef IterativeSolver<SpaceType,  LocalSpaceType> IterativeSolverType;
        typedef CGSolver<SpaceType,  LocalSpaceType> CGSolverType;
        typedef DeflatedCGSolver<SpaceType,  LocalSpaceType> DeflatedCGSolverType;
        typedef BICGSTABSolver<SpaceType,  LocalSpaceType> BICGSTABSolverType;
        typedef TFQMRSolver<SpaceType,  LocalSpaceType> TFQMRSolverType;
        typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType> SkylineLUFactorizationSolverType;
        typedef AMGCLSolver<SpaceType,  LocalSpaceType> AMGCLSolverType;
        typedef AMGCL_NS_Solver<SpaceType,  LocalSpaceType> AMGCL_NS_SolverType;
        typedef SkylineLUCustomScalarSolver<ComplexSpaceType, ComplexLocalSpaceType> SkylineLUComplexSolverType;

        typedef ScalingSolver<SpaceType,  LocalSpaceType> ScalingSolverType;
        typedef MonotonicityPreservingSolver<SpaceType,  LocalSpaceType> MonotonicityPreservingSolverType;

        //NOTE: here we must create persisting objects for the linear solvers
        static auto CGSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,CGSolverType>();
        static auto BICGSTABSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,BICGSTABSolverType>();
        static auto DeflatedCGSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,DeflatedCGSolverType>();
        static auto SkylineLUFactorizationSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,SkylineLUFactorizationSolverType>();
        static auto TFQMRSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,TFQMRSolverType>();
        static auto AMGCLSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,AMGCLSolverType>();
        static auto AMGCL_NS_SolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,AMGCL_NS_SolverType>();
        static auto ScalingSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,ScalingSolverType>();
        static auto MonotonicityPreservingSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,MonotonicityPreservingSolverType>();
        static auto SkylineLUComplexSolverFactory = StandardLinearSolverFactory<ComplexSpaceType, ComplexLocalSpaceType, SkylineLUComplexSolverType>();

        //registration of linear solvers
//         KRATOS_REGISTER_LINEAR_SOLVER("LinearSolver", StandardLinearSolverFactory<SpaceType,LocalSpaceType,LinearSolverType>());
        KRATOS_REGISTER_LINEAR_SOLVER("cg", CGSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("bicgstab", BICGSTABSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("deflated_cg", DeflatedCGSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("tfqmr", TFQMRSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("skyline_lu_factorization",SkylineLUFactorizationSolverFactory );
        KRATOS_REGISTER_LINEAR_SOLVER("amgcl", AMGCLSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("amgcl_ns",AMGCL_NS_SolverFactory );
        KRATOS_REGISTER_LINEAR_SOLVER("scaling",ScalingSolverFactory );
        KRATOS_REGISTER_LINEAR_SOLVER("monotonicity_preserving",MonotonicityPreservingSolverFactory );
        KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("skyline_lu_complex", SkylineLUComplexSolverFactory);

    };
} // Namespace Kratos
