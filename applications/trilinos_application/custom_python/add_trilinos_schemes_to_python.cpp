//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

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
#include "trilinos_application.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

//strategies
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"
#include "custom_strategies/schemes/trilinos_residualbased_newmark_scheme.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_variable_property_static_scheme.h"

// FluidDynamicsApplication schemes
#include "../../FluidDynamicsApplication/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"
#include "../../FluidDynamicsApplication/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bdf_scheme_turbulent.h"
#include "../../FluidDynamicsApplication/custom_strategies/strategies/gear_scheme.h"

// Incompressible fluid
#include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
#include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_lagrangian_monolithic_scheme.h"
#include "../../incompressible_fluid_application/custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme_dpg_enriched.h"

// AdjointFluidApplication
#include "../../AdjointFluidApplication/custom_utilities/response_function.h"
#include "../../AdjointFluidApplication/custom_schemes/adjoint_steady_velocity_pressure_scheme.h"
#include "../../AdjointFluidApplication/custom_schemes/adjoint_bossak_scheme.h"

//teuchos parameter list
#include "Teuchos_ParameterList.hpp"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

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

void  AddSchemes(pybind11::module& m)
{
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef Scheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosBaseSchemeType;
    typedef ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosResidualBasedIncrementalUpdateStaticSchemeType;

//********************************************************************
    //********************************************************************
    class_< TrilinosBaseSchemeType, typename TrilinosBaseSchemeType::Pointer >
    (m, "TrilinosScheme").def(init< >() )
    .def( "Initialize", &TrilinosBaseSchemeType::Initialize )
    .def( "SchemeIsInitialized", &TrilinosBaseSchemeType::SchemeIsInitialized )
    .def( "ElementsAreInitialized", &TrilinosBaseSchemeType::ElementsAreInitialized )
    .def( "ConditionsAreInitialized", &TrilinosBaseSchemeType::ConditionsAreInitialized )
    .def( "InitializeElements", &TrilinosBaseSchemeType::InitializeElements )
    .def( "InitializeConditions", &TrilinosBaseSchemeType::InitializeConditions )
    .def( "InitializeSolutionStep", &TrilinosBaseSchemeType::InitializeSolutionStep )
    .def( "FinalizeSolutionStep", &TrilinosBaseSchemeType::FinalizeSolutionStep )
    .def( "InitializeNonLinIteration", &TrilinosBaseSchemeType::InitializeNonLinIteration )
    .def( "FinalizeNonLinIteration", &TrilinosBaseSchemeType::FinalizeNonLinIteration )
    .def( "Predict", &TrilinosBaseSchemeType::Predict )
    .def( "Update", &TrilinosBaseSchemeType::Update )
    .def( "CalculateOutputData", &TrilinosBaseSchemeType::CalculateOutputData )
    .def( "Clean", &TrilinosBaseSchemeType::Clean )
    .def( "Clear", &TrilinosBaseSchemeType::Clear )
    .def( "MoveMesh", MoveMesh )
    .def("Check", &TrilinosBaseSchemeType::Check )
    ;

    class_ <
        ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >
           (m,"TrilinosResidualBasedIncrementalUpdateStaticScheme")
           .def(init< >() );

    class_ <
        ResidualBasedIncrementalUpdateStaticSchemeSlip< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedIncrementalUpdateStaticSchemeSlip< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>  >
           (
               m,"TrilinosResidualBasedIncrementalUpdateStaticSchemeSlip").def(init< unsigned int, unsigned int >()
           );

    class_ <
        ResidualBasedLagrangianMonolithicScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedLagrangianMonolithicScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >
           (
               m,"TrilinosResidualBasedLagrangianMonolithicScheme").def(init<int >()
           );

    class_ <
        TrilinosResidualBasedNewmarkScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename TrilinosResidualBasedNewmarkScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >
           (
               m,"TrilinosResidualBasedNewmarkScheme").def(init<double >()
           );

    class_ <
        ResidualBasedBossakDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedBossakDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType > (m,"TrilinosResidualBasedBossakDisplacementScheme")
        .def(init<double >())
        ;

    class_ <
        ResidualBasedBDFDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedBDFDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType > (m,"TrilinosResidualBasedBDFDisplacementScheme")
        .def(init<  >())
        .def(init <const std::size_t>())
        ;

    class_ <
        ResidualBasedBDFCustomScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedBDFCustomScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType > (m,"TrilinosResidualBasedBDFCustomScheme")
        .def(init<  >())
        .def(init <const std::size_t>())
        .def(init <const std::size_t, Parameters>())
        ;

    typedef ResidualBasedPredictorCorrectorVelocityBossakScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosResidualBasedPredictorCorrectorVelocityBossak;

    typedef ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TurbulentBossakBaseType;

    class_ <
        TrilinosResidualBasedPredictorCorrectorVelocityBossak,
        typename TrilinosResidualBasedPredictorCorrectorVelocityBossak::Pointer,
        TrilinosBaseSchemeType >
           (
               m,"TrilinosPredictorCorrectorVelocityBossakScheme").def(init<double, double >()
           );

    class_ < TurbulentBossakBaseType, typename TurbulentBossakBaseType::Pointer,TrilinosBaseSchemeType >
        (m,"TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent")
        .def(init<double, double, unsigned int, Process::Pointer >())
        .def(init<double,double,unsigned int >())
        .def(init<double,unsigned int, const Variable<int>&>())
        ;

    class_ <
        ResidualBasedPredictorCorrectorBDFSchemeTurbulent< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedPredictorCorrectorBDFSchemeTurbulent< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >(m,"TrilinosResidualBasedPredictorCorrectorBDFScheme")
        .def(init<unsigned int, Variable<double>& >() );

    class_ <
        ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >
        (m,"TrilinosResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched")
        .def(init<double, double, unsigned int>() );

    class_ <
        TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosResidualBasedIncrementalUpdateStaticSchemeType >
        (m,"TrilinosResidualBasedIncrementalUpdateStaticVariablePropertyScheme")
        .def(init< >());

    typedef GearScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> GearSchemeBaseType;

    class_ < GearSchemeBaseType, typename GearSchemeBaseType::Pointer, TrilinosBaseSchemeType >( m,"TrilinosGearScheme")
            .def(init<Process::Pointer >() )
            .def(init<>()) // constructor without a turbulence model
            .def(init<const Variable<int>&>()) // constructor for periodic conditions
    ;

    class_<
        AdjointSteadyVelocityPressureScheme<TrilinosSparseSpaceType,TrilinosLocalSpaceType>,
        typename AdjointSteadyVelocityPressureScheme<TrilinosSparseSpaceType,TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >( m,"TrilinosAdjointSteadyVelocityPressureScheme")
        .def(init<Parameters&, ResponseFunction::Pointer>() );

    class_<
        AdjointBossakScheme<TrilinosSparseSpaceType,TrilinosLocalSpaceType>,
        typename AdjointBossakScheme<TrilinosSparseSpaceType,TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >( m,"TrilinosAdjointBossakScheme")
        .def(init<Parameters&, ResponseFunction::Pointer>() );
}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
