//
//   Project Name:        Kratos
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-12-09 20:20:55 $
//   Revision:            $Revision: 1.5 $
//
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

#include "custom_python/add_trilinos_strategies_to_python.h"

//Trilinos includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"


// Project includes
#include "includes/define.h"
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/model_part.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
// #include "solving_strategies/schemes/scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
// #include "custom_strategies/schemes/trilinos_residualbased_lagrangian_monolithic_scheme.h"
// #include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
// #include "custom_strategies/schemes/trilinos_predictorcorrector_velocity_bossak_scheme.h"

//convergence criterias
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"
// #include "solving_strategies/convergencecriterias/displacement_criteria.h"
//
// //Builder And Solver
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver_periodic.h"
#include "custom_strategies/builder_and_solvers/trilinos_elimination_builder_and_solver.h"
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"
#include "custom_strategies/strategies/trilinos_convdiff_strategy.h"
#include "custom_strategies/strategies/trilinos_laplacian_meshmoving_strategy.h"
#include "custom_strategies/strategies/trilinos_structural_meshmoving_strategy.h"

//linear solvers
// #include "linear_solvers/linear_solver.h"

//utilities
// #include "python/pointer_vector_set_python_interface.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

// #include "external_includes/aztec_solver.h"
// #include "external_includes/amesos_solver.h"
// #include "external_includes/ml_solver.h"

//configuration files
#include "incompressible_fluid_application/custom_strategies/strategies/solver_configuration.h"
//~ #include "custom_strategies/strategies/trilinos_fractionalstep_configuration.h"
#include "incompressible_fluid_application/custom_strategies/strategies/fractional_step_strategy.h"
#include "incompressible_fluid_application/incompressible_fluid_application.h"
#include "linear_solvers/linear_solver.h"

#include "FluidDynamicsApplication/custom_strategies/strategies/fs_strategy.h"
#include "FluidDynamicsApplication/custom_utilities/solver_settings.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;

void AddStrategies()
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
    typedef ConvergenceCriteria< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosConvergenceCriteria;

    typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;
    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    typedef BuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;

    //********************************************************************
    //********************************************************************
    //Builder and Solver

    class_< TrilinosBuilderAndSolverType::DofsArrayType, boost::noncopyable > ( "DofsArrayType", init<>() );

    // Builder and solver base class
    class_< TrilinosBuilderAndSolverType, boost::noncopyable >
    ( "TrilinosResidualBasedBuilderAndSolver", init<TrilinosLinearSolverType::Pointer> () )
    .def( "SetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::SetCalculateReactionsFlag )
    .def( "GetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::GetCalculateReactionsFlag )
    .def( "SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::SetDofSetIsInitializedFlag )
    .def( "GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::GetDofSetIsInitializedFlag )
    .def( "SetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::SetReshapeMatrixFlag )
    .def( "GetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::GetReshapeMatrixFlag )
    .def( "GetEquationSystemSize", &TrilinosBuilderAndSolverType::GetEquationSystemSize )
    .def( "BuildLHS", &TrilinosBuilderAndSolverType::BuildLHS )
    .def( "BuildRHS", &TrilinosBuilderAndSolverType::BuildRHS )
    .def( "Build", &TrilinosBuilderAndSolverType::Build )
    .def( "SystemSolve", &TrilinosBuilderAndSolverType::SystemSolve )
    .def( "BuildAndSolve", &TrilinosBuilderAndSolverType::BuildAndSolve )
    .def( "BuildRHSAndSolve", &TrilinosBuilderAndSolverType::BuildRHSAndSolve )
    .def( "ApplyDirichletConditions", &TrilinosBuilderAndSolverType::ApplyDirichletConditions )
    .def( "SetUpDofSet", &TrilinosBuilderAndSolverType::SetUpDofSet )
    .def( "GetDofSet", &TrilinosBuilderAndSolverType::GetDofSet, return_internal_reference<>() )
    .def( "SetUpSystem", &TrilinosBuilderAndSolverType::SetUpSystem )
    .def( "ResizeAndInitializeVectors", &TrilinosBuilderAndSolverType::ResizeAndInitializeVectors )
    .def( "InitializeSolutionStep", &TrilinosBuilderAndSolverType::InitializeSolutionStep )
    .def( "FinalizeSolutionStep", &TrilinosBuilderAndSolverType::FinalizeSolutionStep )
    .def( "CalculateReactions", &TrilinosBuilderAndSolverType::CalculateReactions )
    .def( "Clear", &TrilinosBuilderAndSolverType::Clear )
    .def( "SetEchoLevel", &TrilinosBuilderAndSolverType::SetEchoLevel )
    .def( "GetEchoLevel", &TrilinosBuilderAndSolverType::GetEchoLevel )
    ;
    
    typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosResidualBasedEliminationBuilderAndSolverType;
    class_< TrilinosResidualBasedEliminationBuilderAndSolverType, bases<TrilinosBuilderAndSolverType>, boost::noncopyable >
    ( "TrilinosEliminationBuilderAndSolver", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > () )
    ;

    typedef TrilinosBlockBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBlockBuilderAndSolverType;
    class_< TrilinosBlockBuilderAndSolverType,bases< TrilinosBuilderAndSolverType >, boost::noncopyable >
    ( "TrilinosBlockBuilderAndSolver", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > () )
    ;
    
    class_< TrilinosBlockBuilderAndSolverPeriodic< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >,
    bases< TrilinosBlockBuilderAndSolverType >, boost::noncopyable >
    ( "TrilinosBlockBuilderAndSolverPeriodic", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer, Kratos::Variable<int>& >() )
    ;

    //********************************************************************************************
    class_ < SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, boost::noncopyable >
    ( "SolverConfiguration", init< ModelPart&, unsigned int>() )
    .def( "GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag )
    .def( "GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag )
    ;

    //********************************************************************************************
    class_< FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, boost::noncopyable >
    ( "TrilinosFractionalStepStrategy",
      init < ModelPart&,
      SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >&,
      bool,
      double, double,
      int, int,
      unsigned int, unsigned int,
      bool
      > () )
    .def( "Solve", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Solve )
    .def( "SolveStep1", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep1 )
    .def( "SolveStep2", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep2 )
    .def( "SolveStep3", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep3 )
    .def( "SolveStep4", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep4 )
    .def( "ActOnLonelyNodes", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ActOnLonelyNodes )
    .def( "Clear", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Clear )
    .def( "FractionalVelocityIteration", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::FractionalVelocityIteration )
    .def( "ComputeReactions",&FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ComputeReactions)
    .def( "ConvergenceCheck", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ConvergenceCheck )
    .def( "InitializeFractionalStep", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::InitializeFractionalStep )
    .def( "PredictVelocity", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::PredictVelocity )
    .def( "InitializeProjections", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::InitializeProjections )
    .def( "AssignInitialStepValues", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::AssignInitialStepValues )
    .def( "IterativeSolve", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::IterativeSolve )
    .def( "SavePressureIteration", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SavePressureIteration )
    .def( "ApplyFractionalVelocityFixity", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ApplyFractionalVelocityFixity )
    .def( "SetEchoLevel", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SetEchoLevel )
    ;

    //********************************************************************************************
    class_< TrilinosConvectionDiffusionStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, boost::noncopyable >
    ( "TrilinosConvectionDiffusionStrategy", init < Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, bool, int, int > () )
    .def( "Solve", &TrilinosConvectionDiffusionStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Solve )
    ;
    
    // Strategy base class
    class_< TrilinosBaseSolvingStrategyType, boost::noncopyable > ("TrilinosSolvingStrategy", init < ModelPart&, bool >())
    .def("Predict", &TrilinosBaseSolvingStrategyType::Predict)
    .def("Initialize", &TrilinosBaseSolvingStrategyType::Initialize)
    .def("Solve", &TrilinosBaseSolvingStrategyType::Solve)
    .def("IsConverged", &TrilinosBaseSolvingStrategyType::IsConverged)
    .def("CalculateOutputData", &TrilinosBaseSolvingStrategyType::CalculateOutputData)
    .def("SetEchoLevel", &TrilinosBaseSolvingStrategyType::SetEchoLevel)
    .def("GetEchoLevel", &TrilinosBaseSolvingStrategyType::GetEchoLevel)
    .def("SetRebuildLevel", &TrilinosBaseSolvingStrategyType::SetRebuildLevel)
    .def("GetRebuildLevel", &TrilinosBaseSolvingStrategyType::GetRebuildLevel)
    .def("SetMoveMeshFlag", &TrilinosBaseSolvingStrategyType::SetMoveMeshFlag)
    .def("MoveMeshFlag", &TrilinosBaseSolvingStrategyType::MoveMeshFlag)
    .def("MoveMesh", &TrilinosBaseSolvingStrategyType::MoveMesh)
    .def("Clear", &TrilinosBaseSolvingStrategyType::Clear)
    .def("Check", &TrilinosBaseSolvingStrategyType::Check)
    .def("InitializeSolutionStep", &TrilinosBaseSolvingStrategyType::InitializeSolutionStep)
    .def("FinalizeSolutionStep", &TrilinosBaseSolvingStrategyType::FinalizeSolutionStep)
    .def("SolveSolutionStep", &TrilinosBaseSolvingStrategyType::SolveSolutionStep)
    //.def("GetModelPart", &BaseSolvingStrategyType::GetModelPart )
    ;

    typedef FSStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosFSStrategy;
    class_< TrilinosFSStrategy, boost::noncopyable >
    ("TrilinosFSStrategy",init< ModelPart&, SolverSettings< TrilinosSparseSpaceType,TrilinosLocalSpaceType, TrilinosLinearSolverType >&, bool >())
    .def(init< ModelPart&, SolverSettings< TrilinosSparseSpaceType,TrilinosLocalSpaceType, TrilinosLinearSolverType >&, bool, const Kratos::Variable<int>& >())
    .def("Solve",&TrilinosFSStrategy::Solve)
    .def("Clear",&TrilinosFSStrategy::Clear)
    .def("CalculateReactions",&TrilinosFSStrategy::CalculateReactions)
    .def("AddIterationStep",&TrilinosFSStrategy::AddIterationStep)
    .def("ClearExtraIterationSteps",&TrilinosFSStrategy::ClearExtraIterationSteps)
    ;
            
    typedef ResidualBasedNewtonRaphsonStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosNewtonRaphsonStrategy;
    class_< TrilinosNewtonRaphsonStrategy , bases< TrilinosBaseSolvingStrategyType >, boost::noncopyable >
    ("TrilinosNewtonRaphsonStrategy", init< ModelPart&, TrilinosBaseSchemeType::Pointer, TrilinosLinearSolverType::Pointer, TrilinosConvergenceCriteria::Pointer, TrilinosBuilderAndSolverType::Pointer, int, bool, bool, bool >())
    ;
            
    // Mesh Moving ********************************************************************************************
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

    class_< TrilinosLaplacianMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, boost::noncopyable >
    ("TrilinosLaplacianMeshMovingStrategy", init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, int, int, bool >() )
    .def("MoveNodes",&TrilinosLaplacianMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::MoveNodes)
    .def("Solve", &TrilinosLaplacianMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Solve )
    .def("SetEchoLevel", &TrilinosLaplacianMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SetEchoLevel )
    ;

    class_< TrilinosStructuralMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, boost::noncopyable >
    ("TrilinosStructuralMeshMovingStrategy", init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, int, bool >() )
    .def("MoveNodes",&TrilinosStructuralMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::MoveNodes)
    .def("Solve", &TrilinosStructuralMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Solve )
    .def("SetEchoLevel", &TrilinosStructuralMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SetEchoLevel )
    ;

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
