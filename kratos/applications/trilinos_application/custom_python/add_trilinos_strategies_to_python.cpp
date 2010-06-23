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
// #include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver_deactivation.h"
#include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
#include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_vec.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_mixed.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_deactivation.h"
#include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_deactivation_vec.h"
#include "custom_strategies/builder_and_solvers/trilinos_pressure_splitting_builder_and_solver.h"

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
#include "../../incompressible_fluid_application/custom_strategies/strategies/solver_configuration.h"
#include "custom_strategies/strategies/trilinos_fractionalstep_configuration.h"
#include "../../incompressible_fluid_application/custom_strategies/strategies/fractional_step_strategy.h"
#include "../../incompressible_fluid_application/incompressible_fluid_application.h"



namespace Kratos {

    namespace Python {

        using namespace boost::python;


	void  AddStrategies()
        {
       typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
        typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

            typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;
            typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
            typedef TrilinosResidualBasedEliminationBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverType;


            //********************************************************************
            //********************************************************************
            //Builder and Solver

            class_< TrilinosBuilderAndSolverType::DofsArrayType, boost::noncopyable > ("DofsArrayType", init<>());

            class_< TrilinosBuilderAndSolverType, boost::noncopyable >
                    ("TrilinosResidualBasedBuilderAndSolver", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverType::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverType::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverType::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosBuilderAndSolverType::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosBuilderAndSolverType::BuildLHS)
                    .def("BuildRHS", &TrilinosBuilderAndSolverType::BuildRHS)
                    .def("Build", &TrilinosBuilderAndSolverType::Build)
                    .def("SystemSolve", &TrilinosBuilderAndSolverType::SystemSolve)
                    .def("BuildAndSolve", &TrilinosBuilderAndSolverType::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverType::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverType::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosBuilderAndSolverType::SetUpDofSet)
                    .def("GetDofSet", &TrilinosBuilderAndSolverType::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosBuilderAndSolverType::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverType::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosBuilderAndSolverType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverType::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosBuilderAndSolverType::CalculateReactions)
                    .def("Clear", &TrilinosBuilderAndSolverType::Clear)
                    .def("SetEchoLevel", &TrilinosBuilderAndSolverType::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosBuilderAndSolverType::GetEchoLevel)
                    ;

            typedef TrilinosResidualBasedEliminationBuilderAndSolverDeactivation< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverDeactivationType;

            class_< TrilinosBuilderAndSolverDeactivationType, boost::noncopyable >
                    ("TrilinosResidualBasedBuilderAndSolverDeactivation", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverDeactivationType::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverDeactivationType::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverDeactivationType::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverDeactivationType::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverDeactivationType::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverDeactivationType::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosBuilderAndSolverDeactivationType::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosBuilderAndSolverDeactivationType::BuildLHS)
                    .def("BuildRHS", &TrilinosBuilderAndSolverDeactivationType::BuildRHS)
                    .def("Build", &TrilinosBuilderAndSolverDeactivationType::Build)
                    .def("SystemSolve", &TrilinosBuilderAndSolverDeactivationType::SystemSolve)
                    .def("BuildAndSolve", &TrilinosBuilderAndSolverDeactivationType::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverDeactivationType::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverDeactivationType::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosBuilderAndSolverDeactivationType::SetUpDofSet)
                    .def("GetDofSet", &TrilinosBuilderAndSolverDeactivationType::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosBuilderAndSolverDeactivationType::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverDeactivationType::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosBuilderAndSolverDeactivationType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverDeactivationType::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosBuilderAndSolverDeactivationType::CalculateReactions)
                    .def("Clear", &TrilinosBuilderAndSolverDeactivationType::Clear)
                    .def("SetEchoLevel", &TrilinosBuilderAndSolverDeactivationType::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosBuilderAndSolverDeactivationType::GetEchoLevel)
                    ;


            typedef TrilinosBuilderAndSolverML< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverMLtype;

            class_< TrilinosBuilderAndSolverMLtype, bases<TrilinosBuilderAndSolverType>, boost::noncopyable >
                    ("TrilinosBuilderAndSolverML", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLtype::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLtype::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLtype::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLtype::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLtype::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLtype::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosBuilderAndSolverMLtype::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosBuilderAndSolverMLtype::BuildLHS)
                    .def("BuildRHS", &TrilinosBuilderAndSolverMLtype::BuildRHS)
                    .def("Build", &TrilinosBuilderAndSolverMLtype::Build)
                    .def("SystemSolve", &TrilinosBuilderAndSolverMLtype::SystemSolve)
                    .def("BuildAndSolve", &TrilinosBuilderAndSolverMLtype::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverMLtype::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverMLtype::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosBuilderAndSolverMLtype::SetUpDofSet)
                    .def("GetDofSet", &TrilinosBuilderAndSolverMLtype::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosBuilderAndSolverMLtype::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverMLtype::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosBuilderAndSolverMLtype::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverMLtype::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosBuilderAndSolverMLtype::CalculateReactions)
                    .def("Clear", &TrilinosBuilderAndSolverMLtype::Clear)
                    .def("SetEchoLevel", &TrilinosBuilderAndSolverMLtype::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosBuilderAndSolverMLtype::GetEchoLevel)
                    ;

            typedef TrilinosBuilderAndSolverMLDeactivation< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverMLDeactivationtype;

            class_< TrilinosBuilderAndSolverMLDeactivationtype, bases<TrilinosBuilderAndSolverType>, boost::noncopyable >
                    ("TrilinosBuilderAndSolverMLDeactivation", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLDeactivationtype::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLDeactivationtype::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLDeactivationtype::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLDeactivationtype::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLDeactivationtype::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLDeactivationtype::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosBuilderAndSolverMLDeactivationtype::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosBuilderAndSolverMLDeactivationtype::BuildLHS)
                    .def("BuildRHS", &TrilinosBuilderAndSolverMLDeactivationtype::BuildRHS)
                    .def("Build", &TrilinosBuilderAndSolverMLDeactivationtype::Build)
                    .def("SystemSolve", &TrilinosBuilderAndSolverMLDeactivationtype::SystemSolve)
                    .def("BuildAndSolve", &TrilinosBuilderAndSolverMLDeactivationtype::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverMLDeactivationtype::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverMLDeactivationtype::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosBuilderAndSolverMLDeactivationtype::SetUpDofSet)
                    .def("GetDofSet", &TrilinosBuilderAndSolverMLDeactivationtype::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosBuilderAndSolverMLDeactivationtype::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverMLDeactivationtype::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosBuilderAndSolverMLDeactivationtype::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverMLDeactivationtype::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosBuilderAndSolverMLDeactivationtype::CalculateReactions)
                    .def("Clear", &TrilinosBuilderAndSolverMLDeactivationtype::Clear)
                    .def("SetEchoLevel", &TrilinosBuilderAndSolverMLDeactivationtype::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosBuilderAndSolverMLDeactivationtype::GetEchoLevel)
                    ;



            typedef TrilinosBuilderAndSolverML2D< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverML2Dtype;

            class_< TrilinosBuilderAndSolverML2Dtype, boost::noncopyable >
                    ("TrilinosBuilderAndSolverML2D", init<Epetra_MpiComm&, int, int, TrilinosLinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverML2Dtype::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverML2Dtype::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverML2Dtype::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverML2Dtype::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverML2Dtype::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverML2Dtype::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosBuilderAndSolverML2Dtype::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosBuilderAndSolverML2Dtype::BuildLHS)
                    .def("BuildRHS", &TrilinosBuilderAndSolverML2Dtype::BuildRHS)
                    .def("Build", &TrilinosBuilderAndSolverML2Dtype::Build)
                    .def("SystemSolve", &TrilinosBuilderAndSolverML2Dtype::SystemSolve)
                    .def("BuildAndSolve", &TrilinosBuilderAndSolverML2Dtype::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverML2Dtype::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverML2Dtype::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosBuilderAndSolverML2Dtype::SetUpDofSet)
                    .def("GetDofSet", &TrilinosBuilderAndSolverML2Dtype::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosBuilderAndSolverML2Dtype::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverML2Dtype::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosBuilderAndSolverML2Dtype::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverML2Dtype::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosBuilderAndSolverML2Dtype::CalculateReactions)
                    .def("Clear", &TrilinosBuilderAndSolverML2Dtype::Clear)
                    .def("SetEchoLevel", &TrilinosBuilderAndSolverML2Dtype::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosBuilderAndSolverML2Dtype::GetEchoLevel)
                    ;


            typedef TrilinosBuilderAndSolverMLmixed< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverMLmixedType;

            class_< TrilinosBuilderAndSolverMLmixedType, boost::noncopyable >
                    ("TrilinosBuilderAndSolverMLmixed", init<Epetra_MpiComm&, int, int, TrilinosLinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLmixedType::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLmixedType::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLmixedType::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLmixedType::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLmixedType::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLmixedType::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosBuilderAndSolverMLmixedType::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosBuilderAndSolverMLmixedType::BuildLHS)
                    .def("BuildRHS", &TrilinosBuilderAndSolverMLmixedType::BuildRHS)
                    .def("Build", &TrilinosBuilderAndSolverMLmixedType::Build)
                    .def("SystemSolve", &TrilinosBuilderAndSolverMLmixedType::SystemSolve)
                    .def("BuildAndSolve", &TrilinosBuilderAndSolverMLmixedType::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverMLmixedType::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverMLmixedType::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosBuilderAndSolverMLmixedType::SetUpDofSet)
                    .def("GetDofSet", &TrilinosBuilderAndSolverMLmixedType::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosBuilderAndSolverMLmixedType::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverMLmixedType::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosBuilderAndSolverMLmixedType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverMLmixedType::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosBuilderAndSolverMLmixedType::CalculateReactions)
                    .def("Clear", &TrilinosBuilderAndSolverMLmixedType::Clear)
                    .def("SetEchoLevel", &TrilinosBuilderAndSolverMLmixedType::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosBuilderAndSolverMLmixedType::GetEchoLevel)
                    ;

            typedef TrilinosBuilderAndSolverMLDeactivation2D< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBuilderAndSolverMLDeactivation2Dtype;

            class_< TrilinosBuilderAndSolverMLDeactivation2Dtype, boost::noncopyable >
                    ("TrilinosBuilderAndSolverMLDeactivation2D", init<Epetra_MpiComm&, int, int, TrilinosLinearSolverType::Pointer > ())
                    .def("SetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLDeactivation2Dtype::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosBuilderAndSolverMLDeactivation2Dtype::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLDeactivation2Dtype::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosBuilderAndSolverMLDeactivation2Dtype::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLDeactivation2Dtype::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosBuilderAndSolverMLDeactivation2Dtype::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosBuilderAndSolverMLDeactivation2Dtype::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosBuilderAndSolverMLDeactivation2Dtype::BuildLHS)
                    .def("BuildRHS", &TrilinosBuilderAndSolverMLDeactivation2Dtype::BuildRHS)
                    .def("Build", &TrilinosBuilderAndSolverMLDeactivation2Dtype::Build)
                    .def("SystemSolve", &TrilinosBuilderAndSolverMLDeactivation2Dtype::SystemSolve)
                    .def("BuildAndSolve", &TrilinosBuilderAndSolverMLDeactivation2Dtype::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosBuilderAndSolverMLDeactivation2Dtype::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosBuilderAndSolverMLDeactivation2Dtype::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosBuilderAndSolverMLDeactivation2Dtype::SetUpDofSet)
                    .def("GetDofSet", &TrilinosBuilderAndSolverMLDeactivation2Dtype::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosBuilderAndSolverMLDeactivation2Dtype::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosBuilderAndSolverMLDeactivation2Dtype::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosBuilderAndSolverMLDeactivation2Dtype::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosBuilderAndSolverMLDeactivation2Dtype::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosBuilderAndSolverMLDeactivation2Dtype::CalculateReactions)
                    .def("Clear", &TrilinosBuilderAndSolverMLDeactivation2Dtype::Clear)
                    .def("SetEchoLevel", &TrilinosBuilderAndSolverMLDeactivation2Dtype::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosBuilderAndSolverMLDeactivation2Dtype::GetEchoLevel)
                    ;

            typedef TrilinosPressureSplittingBuilderAndSolver< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >
                    TrilinosPressureSplittingBuilderAndSolverType;

            class_< TrilinosPressureSplittingBuilderAndSolverType, boost::noncopyable >
                    ("TrilinosPressureSplittingBuilderAndSolver", init<Epetra_MpiComm&, int, TrilinosLinearSolverType::Pointer, TrilinosLinearSolverType::Pointer, unsigned int, bool, double, double, double > ())
                    .def("SetCalculateReactionsFlag", &TrilinosPressureSplittingBuilderAndSolverType::SetCalculateReactionsFlag)
                    .def("GetCalculateReactionsFlag", &TrilinosPressureSplittingBuilderAndSolverType::GetCalculateReactionsFlag)
                    .def("SetDofSetIsInitializedFlag", &TrilinosPressureSplittingBuilderAndSolverType::SetDofSetIsInitializedFlag)
                    .def("GetDofSetIsInitializedFlag", &TrilinosPressureSplittingBuilderAndSolverType::GetDofSetIsInitializedFlag)
                    .def("SetReshapeMatrixFlag", &TrilinosPressureSplittingBuilderAndSolverType::SetReshapeMatrixFlag)
                    .def("GetReshapeMatrixFlag", &TrilinosPressureSplittingBuilderAndSolverType::GetReshapeMatrixFlag)
                    .def("GetEquationSystemSize", &TrilinosPressureSplittingBuilderAndSolverType::GetEquationSystemSize)
                    .def("BuildLHS", &TrilinosPressureSplittingBuilderAndSolverType::BuildLHS)
                    .def("BuildRHS", &TrilinosPressureSplittingBuilderAndSolverType::BuildRHS)
                    .def("Build", &TrilinosPressureSplittingBuilderAndSolverType::Build)
                    .def("SystemSolve", &TrilinosPressureSplittingBuilderAndSolverType::SystemSolve)
                    .def("BuildAndSolve", &TrilinosPressureSplittingBuilderAndSolverType::BuildAndSolve)
                    .def("BuildRHSAndSolve", &TrilinosPressureSplittingBuilderAndSolverType::BuildRHSAndSolve)
                    .def("ApplyDirichletConditions", &TrilinosPressureSplittingBuilderAndSolverType::ApplyDirichletConditions)
                    .def("SetUpDofSet", &TrilinosPressureSplittingBuilderAndSolverType::SetUpDofSet)
                    .def("GetDofSet", &TrilinosPressureSplittingBuilderAndSolverType::GetDofSet, return_internal_reference<>())
                    .def("SetUpSystem", &TrilinosPressureSplittingBuilderAndSolverType::SetUpSystem)
                    .def("ResizeAndInitializeVectors", &TrilinosPressureSplittingBuilderAndSolverType::ResizeAndInitializeVectors)
                    .def("InitializeSolutionStep", &TrilinosPressureSplittingBuilderAndSolverType::InitializeSolutionStep)
                    .def("FinalizeSolutionStep", &TrilinosPressureSplittingBuilderAndSolverType::FinalizeSolutionStep)
                    .def("CalculateReactions", &TrilinosPressureSplittingBuilderAndSolverType::CalculateReactions)
                    .def("Clear", &TrilinosPressureSplittingBuilderAndSolverType::Clear)
                    .def("SetEchoLevel", &TrilinosPressureSplittingBuilderAndSolverType::SetEchoLevel)
                    .def("GetEchoLevel", &TrilinosPressureSplittingBuilderAndSolverType::GetEchoLevel)
                    ;


            //********************************************************************************************
            class_< SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >,
                    boost::noncopyable >
                    ("SolverConfiguration", init< ModelPart&, unsigned int>())
                    .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
                    .def("GetActualizeRHSflag", &ConvergenceCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType >::GetActualizeRHSflag)
                    ;

            class_< TrilinosFractionalStepConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >,
                    bases< SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > >,
                    boost::noncopyable >
                    ("TrilinosFractionalStepConfiguration", init < Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, TrilinosLinearSolverType::Pointer,
                    unsigned int, unsigned int, bool >());


            //********************************************************************************************
            class_< FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >, boost::noncopyable >
                    ("TrilinosFractionalStepStrategy",
                    init < ModelPart&,
                    SolverConfiguration<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >&,
                    bool,
                    double, double,
                    int, int,
                    unsigned int, unsigned int,
                    bool
                    >())
                    .def("Solve", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Solve)
                    .def("SolveStep1", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep1)
                    .def("SolveStep2", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep2)
                    .def("SolveStep3", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep3)
                    .def("SolveStep4", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SolveStep4)
                    .def("ActOnLonelyNodes", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ActOnLonelyNodes)
                    .def("Clear", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::Clear)
                    .def("FractionalVelocityIteration", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::FractionalVelocityIteration)
                    .def("ConvergenceCheck", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ConvergenceCheck)
                    .def("InitializeFractionalStep", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::InitializeFractionalStep)
                    .def("PredictVelocity", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::PredictVelocity)
                    .def("InitializeProjections", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::InitializeProjections)
                    .def("AssignInitialStepValues", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::AssignInitialStepValues)
                    .def("IterativeSolve", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::IterativeSolve)
                    .def("SavePressureIteration", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SavePressureIteration)
                    .def("ApplyFractionalVelocityFixity", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::ApplyFractionalVelocityFixity)
                    .def("SetEchoLevel", &FractionalStepStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >::SetEchoLevel )
            ;

        }


    } // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
