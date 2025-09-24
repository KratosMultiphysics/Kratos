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

#include "custom_elements/membrane_elements/membrane_element.hpp"
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
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    Element                                             element;

    Matrix                   LHS_Contribution;
    Vector                   RHS_Contribution;
    std::vector<std::size_t> EquationId;
    ProcessInfo              CurrentProcessInfo;

    scheme.CalculateRHSContribution(element, RHS_Contribution, EquationId, CurrentProcessInfo);
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeCanInitializeSolutionStep, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    Element                                             element;

    Matrix                   LHS_Contribution;
    Vector                   RHS_Contribution;
    std::vector<std::size_t> EquationId;
    ProcessInfo              CurrentProcessInfo;

    Model            model;
    auto&            model_part = model.CreateModelPart("Main");
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    scheme.InitializeSolutionStep(model_part, A, Dx, b);
}

class MockElementForLoadSteppingScheme : public Element
{
public:
    void Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == INTERNAL_FORCES_VECTOR) Output = ScalarVector(10, 1);
    }
};

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSIsEqualToCurrentInternalForcesAtStartOfStage,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    MockElementForLoadSteppingScheme                    element;

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 0.0;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model            model;
    auto&            model_part = model.CreateModelPart("Main");
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;
    scheme.CalculateRHSContribution(element, actual_right_hand_side, EquationId, CurrentProcessInfo);

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, Vector{ScalarVector(10, 1)}, 1e-6);
}

} // namespace Kratos::Testing