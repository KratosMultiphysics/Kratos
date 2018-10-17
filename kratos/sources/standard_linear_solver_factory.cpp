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

        typedef PreconditionerFactoryBase<SpaceType,  LocalSpaceType> PreconditionerFactoryBaseType;

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
        KratosComponents<PreconditionerFactoryBaseType>::Add("None", PreconditionerF);
        KratosComponents<PreconditionerFactoryBaseType>::Add("DiagonalPreconditioner", DiagonalPreconditionerFactory);
        KratosComponents<PreconditionerFactoryBaseType>::Add("ILU0Preconditioner", ILU0Factory);
        KratosComponents<PreconditionerFactoryBaseType>::Add("ILUPreconditioner",ILUFactory );

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

        typedef LinearSolverFactoryBase<SpaceType,  LocalSpaceType> LinearSolverFactoryBaseType;
        typedef LinearSolverFactoryBase<ComplexSpaceType,  ComplexLocalSpaceType> ComplexLinearSolverFactoryBaseType;

        //NOTE: here we must create persisting objects for the linear solvers
        static auto CGSolverFactory = LinearSolverFactory<SpaceType,LocalSpaceType,CGSolverType>();
        static auto BICGSTABSolverFactory = LinearSolverFactory<SpaceType,LocalSpaceType,BICGSTABSolverType>();
        static auto DeflatedCGSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,DeflatedCGSolverType>();
        static auto SkylineLUFactorizationSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,SkylineLUFactorizationSolverType>();
//         static auto MixedUPLinearSolverFactory = LinearSolverFactory<SpaceType,LocalSpaceType,MixedUPLinearSolverType>();
        static auto TFQMRSolverFactory = LinearSolverFactory<SpaceType,LocalSpaceType,TFQMRSolverType>();
        static auto AMGCLSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,AMGCLSolverType>();
        static auto AMGCL_NS_SolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,AMGCL_NS_SolverType>();
        static auto ScalingSolverFactory= LinearSolverFactory<SpaceType,LocalSpaceType,ScalingSolverType>();
        static auto SkylineLUComplexSolverFactory = LinearSolverFactory<ComplexSpaceType, ComplexLocalSpaceType, SkylineLUComplexSolverType>();

        //registration of linear solvers
//         KratosComponents<LinearSolverFactoryBaseType>::Add("LinearSolver", LinearSolverFactory<SpaceType,LocalSpaceType,LinearSolverType>());
        KratosComponents<LinearSolverFactoryBaseType>::Add("CGSolver", CGSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add("BICGSTABSolver", BICGSTABSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add("DeflatedCGSolver", DeflatedCGSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add("TFQMRSolver", TFQMRSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add("SkylineLUFactorizationSolver",SkylineLUFactorizationSolverFactory );
        KratosComponents<LinearSolverFactoryBaseType>::Add("AMGCL", AMGCLSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add("AMGCLSolver", AMGCLSolverFactory); //registered with two different names
        KratosComponents<LinearSolverFactoryBaseType>::Add("AMGCL_NS_Solver",AMGCL_NS_SolverFactory );
        KratosComponents<LinearSolverFactoryBaseType>::Add("ScalingSolver",ScalingSolverFactory );
        KratosComponents<ComplexLinearSolverFactoryBaseType>::Add("SkylineLUComplexSolver", SkylineLUComplexSolverFactory);
        KratosComponents<ComplexLinearSolverFactoryBaseType>::Add("complex_skyline_lu_solver", SkylineLUComplexSolverFactory); // NOTE: Name duplicated for retrocompatibility

    };




} // Namespace Kratos

