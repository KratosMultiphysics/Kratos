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

#include "custom_python/add_trilinos_schemes_to_python.h"

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
// #include "add_trilinos_linear_solvers_to_python.h"
#include "includes/model_part.h"

//strategies
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/schemes/trilinos_residualbased_lagrangian_monolithic_scheme.h"
#include "custom_strategies/schemes/trilinos_residualbased_newmark_scheme.h"
#include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
#include "custom_strategies/schemes/trilinos_predictorcorrector_velocity_bossak_scheme.h"
#include "../../FluidDynamicsApplication/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "custom_strategies/schemes/trilinos_predictorcorrector_velocity_bossak_scheme_turbulent.h"
//#include "../../FluidDynamicsApplication/custom_strategies/strategies/gear_scheme.h"
//#include "custom_strategies/schemes/trilinos_gear_scheme.h"

//convergence criterias
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"
// #include "solving_strategies/convergencecriterias/displacement_criteria.h"
//
// //Builder And Solver
// // #include "solving_strategies/builder_and_solvers/builder_and_solver.h"
// #include "custom_strategies/builder_and_solvers/trilinos_residualbased_elimination_builder_and_solver.h"
// #include "custom_strategies/convergencecriterias/trilinos_displacement_criteria.h"
// #include "custom_strategies/convergencecriterias/trilinos_up_criteria.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_vec.h"
// #include "custom_strategies/builder_and_solvers/trilinos_builder_and_solver_ML_mixed.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//utilities
#include "python/pointer_vector_set_python_interface.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"

#include "external_includes/aztec_solver.h"
#include "external_includes/amesos_solver.h"
#include "external_includes/ml_solver.h"

//configuration files
// #include "../../incompressible_fluid_application/custom_strategies/strategies/solver_configuration.h"
// #include "custom_strategies/strategies/trilinos_fractionalstep_configuration.h"
// #include "../../incompressible_fluid_application/custom_strategies/strategies/fractional_step_strategy.h"
// #include "../../incompressible_fluid_application/incompressible_fluid_application.h"



namespace Kratos
{

    namespace Python
    {

        using namespace boost::python;

        typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

        void MoveMesh( Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType >& dummy, ModelPart::NodesContainerType& rNodes )
        {
            for ( ModelPart::NodeIterator i = rNodes.begin(); i != rNodes.end(); ++i )
            {
                const array_1d<double, 3 > & disp = i->FastGetSolutionStepValue( DISPLACEMENT );
                ( i )->X() = ( i )->X0() + disp[0];
                ( i )->Y() = ( i )->Y0() + disp[1];
                ( i )->Z() = ( i )->Z0() + disp[2];
            }
        }

        void  AddSchemes()
        {
            typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
            typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
            typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;

//********************************************************************
            //********************************************************************
            class_< TrilinosBaseSchemeType, boost::noncopyable >
            ( "TrilinosScheme", init< >() )
            .def( "Initialize", &TrilinosBaseSchemeType::Initialize )
            .def( "SchemeIsInitialized", &TrilinosBaseSchemeType::SchemeIsInitialized )
            .def( "ElementsAreInitialized", &TrilinosBaseSchemeType::ElementsAreInitialized )
            .def( "InitializeElements", &TrilinosBaseSchemeType::InitializeElements )
            .def( "InitializeSolutionStep", &TrilinosBaseSchemeType::InitializeSolutionStep )
            .def( "FinalizeSolutionStep", &TrilinosBaseSchemeType::FinalizeSolutionStep )
            .def( "InitializeNonLinIteration", &TrilinosBaseSchemeType::InitializeNonLinIteration )
            .def( "FinalizeNonLinIteration", &TrilinosBaseSchemeType::FinalizeNonLinIteration )
            .def( "Predict", &TrilinosBaseSchemeType::Predict )
            .def( "Update", &TrilinosBaseSchemeType::Update )
            .def( "CalculateOutputData", &TrilinosBaseSchemeType::CalculateOutputData )
            .def( "Clean", &TrilinosBaseSchemeType::Clean )
            .def( "MoveMesh", MoveMesh )
            .def("Check", &TrilinosBaseSchemeType::Check )
            ;

            class_ < TrilinosResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
            bases< TrilinosBaseSchemeType >, boost::noncopyable >
            (
                "TrilinosResidualBasedIncrementalUpdateStaticScheme", init< >()
            );

            class_ < TrilinosResidualBasedLagrangianMonolithicScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
            bases< TrilinosBaseSchemeType >, boost::noncopyable >
            (
                "TrilinosResidualBasedLagrangianMonolithicScheme", init<int >()
            );

            class_ < TrilinosResidualBasedNewmarkScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
            bases< TrilinosBaseSchemeType >, boost::noncopyable >
            (
                "TrilinosResidualBasedNewmarkScheme", init<double >()
            );

            typedef ResidualBasedPredictorCorrectorVelocityBossakScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme;

            class_ < TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme,
            bases< TrilinosBaseSchemeType >, boost::noncopyable >
            (
                "TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme", init<double, double >()
            );

            class_ < TrilinosPredictorCorrectorVelocityBossakScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
            bases< TrilinosResidualBasedPredictorCorrectorVelocityBossak_BaseScheme >, boost::noncopyable >
            (
                "TrilinosPredictorCorrectorVelocityBossakScheme", init<double, double >()
            );

            typedef ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TurbulentBossakBaseType;

            class_ < TurbulentBossakBaseType,
            bases< TrilinosBaseSchemeType >, boost::noncopyable >
            (
                "TurbulentBossakBaseType", init<double, double, unsigned int, Process::Pointer >()
            );

            class_ < TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
            bases< TurbulentBossakBaseType >, boost::noncopyable >
            (
                "TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent", init<double, double, unsigned int, Process::Pointer >()
            )
                    .def(init<double,double,unsigned int >())// constructor without a turbulence model
                    ;

//            typedef GearScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> GearSchemeBaseType;
//
//            class_ < GearSchemeBaseType,
//            bases< TrilinosBaseSchemeType >, boost::noncopyable >
//            (
//                "GearSchemeBaseType", init<Process::Pointer >()
//            );
//
//            class_< TrilinosGearScheme<TrilinosSparseSpaceType,TrilinosLocalSpaceType>,
//                    bases<GearSchemeBaseType>, boost::noncopyable >
//            ( "TrilinosGearScheme", init<Process::Pointer>() )
//                    .def(init<>()) // constructor without a turbulence model
//                    ;
        }


    } // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
