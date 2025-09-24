// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_strategies/schemes/load_stepping_scheme.hpp"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"



namespace Kratos::Testing
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;


KRATOS_TEST_CASE_IN_SUITE(CanCreateLoadSteppingScheme, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeCanCalculateRHSContribution, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme <SparseSpaceType, LocalSpaceType> scheme;
    Element element;

    Matrix LHS_Contribution;
    Vector RHS_Contribution;
    std::vector<std::size_t> EquationId;
    ProcessInfo CurrentProcessInfo;

    scheme.CalculateRHSContribution(element, RHS_Contribution, EquationId, CurrentProcessInfo);
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeCanInitializeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme <SparseSpaceType, LocalSpaceType> scheme;
    Element element;

    Matrix LHS_Contribution;
    Vector RHS_Contribution;
    std::vector<std::size_t> EquationId;
    ProcessInfo CurrentProcessInfo;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    CompressedMatrix         A;
    Vector                   Dx;
    Vector                   b;

    scheme.InitializeSolutionStep(model_part, A, Dx, b);
}




} // namespace Kratos::Testing