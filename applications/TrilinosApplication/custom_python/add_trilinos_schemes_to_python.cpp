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

// Project includes
#include "includes/define_python.h"

#include "custom_python/add_trilinos_schemes_to_python.h"

// Trilinos includes
#include "Epetra_FEVector.h"

// Project includes
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

// Schemes
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_steady_scheme.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"

// Response function
#include "response_functions/adjoint_response_function.h"

// //teuchos parameter list
// #include "Teuchos_ParameterList.hpp"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

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

//********************************************************************
    //********************************************************************
    py::class_< TrilinosBaseSchemeType, typename TrilinosBaseSchemeType::Pointer >
    (m, "TrilinosScheme").def(py::init< >() )
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

    py::class_ <
        ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType >
           (m,"TrilinosResidualBasedIncrementalUpdateStaticScheme")
           .def(py::init< >() );

    py::class_ <
        ResidualBasedIncrementalUpdateStaticSchemeSlip< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedIncrementalUpdateStaticSchemeSlip< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        ResidualBasedIncrementalUpdateStaticScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>  >
           (
               m,"TrilinosResidualBasedIncrementalUpdateStaticSchemeSlip").def(py::init< unsigned int, unsigned int >()
           );

    py::class_ <
        ResidualBasedBossakDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedBossakDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType > (m,"TrilinosResidualBasedBossakDisplacementScheme")
        .def(py::init<double >())
        ;

    py::class_ <
        ResidualBasedBDFDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedBDFDisplacementScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType > (m,"TrilinosResidualBasedBDFDisplacementScheme")
        .def(py::init<  >())
        .def(py::init <const std::size_t>())
        ;

    py::class_ <
        ResidualBasedBDFCustomScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>,
        typename ResidualBasedBDFCustomScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType>::Pointer,
        TrilinosBaseSchemeType > (m,"TrilinosResidualBasedBDFCustomScheme")
        .def(py::init<  >())
        .def(py::init <const std::size_t>())
        .def(py::init <const std::size_t, Parameters>())
        ;

    typedef ResidualBasedAdjointStaticScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosResidualBasedAdjointStaticSchemeType;
    py::class_<TrilinosResidualBasedAdjointStaticSchemeType, typename TrilinosResidualBasedAdjointStaticSchemeType::Pointer, TrilinosBaseSchemeType>
        (m, "TrilinosResidualBasedAdjointStaticScheme")
        .def(py::init<AdjointResponseFunction::Pointer>())
    ;

    typedef ResidualBasedAdjointSteadyScheme<TrilinosSparseSpaceType, TrilinosLocalSpaceType> TrilinosResidualBasedAdjointSteadySchemeType;
    py::class_<TrilinosResidualBasedAdjointSteadySchemeType, typename TrilinosResidualBasedAdjointSteadySchemeType::Pointer, TrilinosResidualBasedAdjointStaticSchemeType>
        (m, "TrilinosResidualBasedAdjointSteadyScheme")
        .def(py::init<AdjointResponseFunction::Pointer>())
    ;

    typedef ResidualBasedAdjointBossakScheme< TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosResidualBasedAdjointBossakSchemeType;
    py::class_<TrilinosResidualBasedAdjointBossakSchemeType, typename TrilinosResidualBasedAdjointBossakSchemeType::Pointer, TrilinosBaseSchemeType>
        (m, "TrilinosResidualBasedAdjointBossakScheme")
        .def(py::init<Kratos::Parameters, AdjointResponseFunction::Pointer>())
    ;

}

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON defined
