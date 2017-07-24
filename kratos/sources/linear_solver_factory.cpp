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
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/cg_solver.h"
#include "linear_solvers/deflated_cg_solver.h"
#include "linear_solvers/bicgstab_solver.h"
#include "linear_solvers/tfqmr_solver.h"
#include "linear_solvers/mixedup_linear_solver.h"
#include "linear_solvers/amgcl_solver.h"
#include "linear_solvers/amgcl_ns_solver.h"
#include "linear_solvers/scaling_solver.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"

namespace Kratos
{
    void RegisterPreconditioners()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        
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
        KratosComponents<PreconditionerFactoryBaseType>::Add(std::string("None"), PreconditionerF);
        KratosComponents<PreconditionerFactoryBaseType>::Add(std::string("DiagonalPreconditioner"), DiagonalPreconditionerFactory);
        KratosComponents<PreconditionerFactoryBaseType>::Add(std::string("ILU0Preconditioner"), ILU0Factory);
        KratosComponents<PreconditionerFactoryBaseType>::Add(std::string("ILUPreconditioner"),ILUFactory );

    };    
    
    void RegisterLinearSolvers()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
        typedef IterativeSolver<SpaceType,  LocalSpaceType> IterativeSolverType;
        typedef CGSolver<SpaceType,  LocalSpaceType> CGSolverType;
        typedef DeflatedCGSolver<SpaceType,  LocalSpaceType> DeflatedCGSolverType;
//         typedef MixedUPLinearSolver<SpaceType,  LocalSpaceType> MixedUPLinearSolverType;
        typedef BICGSTABSolver<SpaceType,  LocalSpaceType> BICGSTABSolverType;
        typedef TFQMRSolver<SpaceType,  LocalSpaceType> TFQMRSolverType;
        typedef SkylineLUFactorizationSolver<SpaceType,  LocalSpaceType> SkylineLUFactorizationSolverType;
        typedef AMGCLSolver<SpaceType,  LocalSpaceType> AMGCLSolverType;
        typedef AMGCL_NS_Solver<SpaceType,  LocalSpaceType> AMGCL_NS_SolverType;
        
        
        typedef ScalingSolver<SpaceType,  LocalSpaceType> ScalingSolverType;
       
        typedef LinearSolverFactoryBase<SpaceType,  LocalSpaceType> LinearSolverFactoryBaseType;
        
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
        
        //registration of linear solvers
//         KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("LinearSolver"), LinearSolverFactory<SpaceType,LocalSpaceType,LinearSolverType>());
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("CGSolver"), CGSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("BICGSTABSolver"), BICGSTABSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("DeflatedCGSolver"), DeflatedCGSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("TFQMRSolver"), TFQMRSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("SkylineLUFactorizationSolver"),SkylineLUFactorizationSolverFactory );
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("AMGCL"), AMGCLSolverFactory);
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("AMGCLSolver"), AMGCLSolverFactory); //registered with two different names
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("AMGCL_NS_Solver"),AMGCL_NS_SolverFactory );
        KratosComponents<LinearSolverFactoryBaseType>::Add(std::string("ScalingSolver"),ScalingSolverFactory );

    };
    
    


} // Namespace Kratos

