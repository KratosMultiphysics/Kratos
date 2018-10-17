//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/standard_linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "linear_solvers/mixedup_linear_solver.h"
#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/amgcl_ns_solver.h"
#include "linear_solvers/scaling_solver.h"
#include "linear_solvers/skyline_lu_custom_scalar_solver.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"

namespace Kratos
{
    void RegisterPreconditioners()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;

        typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
        typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
        typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
        typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;

        //NOTE: here we must create persisting objects for the linear solvers
        static auto PreconditionerF = PreconditionerFactory<SpaceType,LocalSpaceType,PreconditionerType>();
        static auto DiagonalPreconditionerFactory = PreconditionerFactory<SpaceType,LocalSpaceType,DiagonalPreconditionerType>();
        static auto ILU0Factory= PreconditionerFactory<SpaceType,LocalSpaceType,ILU0PreconditionerType>();
        static auto ILUFactory= PreconditionerFactory<SpaceType,LocalSpaceType,ILUPreconditionerType>();

        //registration of linear solvers
        KRATOS_REGISTER_PRECONDITIONER("None", PreconditionerF);
        KRATOS_REGISTER_PRECONDITIONER("DiagonalPreconditioner", DiagonalPreconditionerFactory);
        KRATOS_REGISTER_PRECONDITIONER("ILU0Preconditioner", ILU0Factory);
        KRATOS_REGISTER_PRECONDITIONER("ILUPreconditioner",ILUFactory );
    };

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
//         typedef MixedUPLinearSolver<SpaceType,  LocalSpaceType> MixedUPLinearSolverType;
        typedef BICGSTABSolver<SpaceType,  LocalSpaceType> BICGSTABSolverType;
        typedef TFQMRSolver<SpaceType,  LocalSpaceType> TFQMRSolverType;
        typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType> SkylineLUFactorizationSolverType;
        typedef AMGCLSolver<SpaceType,  LocalSpaceType> AMGCLSolverType;
        typedef AMGCL_NS_Solver<SpaceType,  LocalSpaceType> AMGCL_NS_SolverType;
        typedef SkylineLUCustomScalarSolver<ComplexSpaceType, ComplexLocalSpaceType> SkylineLUComplexSolverType;

        typedef ScalingSolver<SpaceType,  LocalSpaceType> ScalingSolverType;

        //NOTE: here we must create persisting objects for the linear solvers
        static auto CGSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,CGSolverType>();
        static auto BICGSTABSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,BICGSTABSolverType>();
        static auto DeflatedCGSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,DeflatedCGSolverType>();
        static auto SkylineLUFactorizationSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,SkylineLUFactorizationSolverType>();
//         static auto MixedUPStandardLinearSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,MixedUPLinearSolverType>();
        static auto TFQMRSolverFactory = StandardLinearSolverFactory<SpaceType,LocalSpaceType,TFQMRSolverType>();
        static auto AMGCLSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,AMGCLSolverType>();
        static auto AMGCL_NS_SolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,AMGCL_NS_SolverType>();
        static auto ScalingSolverFactory= StandardLinearSolverFactory<SpaceType,LocalSpaceType,ScalingSolverType>();
        static auto SkylineLUComplexSolverFactory = StandardLinearSolverFactory<ComplexSpaceType, ComplexLocalSpaceType, SkylineLUComplexSolverType>();

        //registration of linear solvers
//         KRATOS_REGISTER_LINEAR_SOLVER("LinearSolver", StandardLinearSolverFactory<SpaceType,LocalSpaceType,LinearSolverType>());
        KRATOS_REGISTER_LINEAR_SOLVER("CGSolver", CGSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("BICGSTABSolver", BICGSTABSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("DeflatedCGSolver", DeflatedCGSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("TFQMRSolver", TFQMRSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("SkylineLUFactorizationSolver",SkylineLUFactorizationSolverFactory );
        KRATOS_REGISTER_LINEAR_SOLVER("AMGCL", AMGCLSolverFactory);
        KRATOS_REGISTER_LINEAR_SOLVER("AMGCLSolver", AMGCLSolverFactory); //registered with two different names
        KRATOS_REGISTER_LINEAR_SOLVER("AMGCL_NS_Solver",AMGCL_NS_SolverFactory );
        KRATOS_REGISTER_LINEAR_SOLVER("ScalingSolver",ScalingSolverFactory );
        KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("SkylineLUComplexSolver", SkylineLUComplexSolverFactory);
        KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER("complex_skyline_lu_solver", SkylineLUComplexSolverFactory); // NOTE: Name duplicated for retrocompatibility

    };




} // Namespace Kratos

