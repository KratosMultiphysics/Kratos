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

#include "boost/numeric/ublas/assignment.hpp"
#include "custom_strategies/schemes/load_stepping_scheme.hpp"
#include "includes/condition.h"
#include "includes/element.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

class MockElementForLoadSteppingScheme : public Element
{
public:
    void Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == INTERNAL_FORCES_VECTOR) Output = mInternalForces;
        if (rVariable == EXTERNAL_FORCES_VECTOR) Output = mExternalForces;
    }

    void SetInternalForces(const Vector& rInternalForces) { mInternalForces = rInternalForces; }

    void SetExternalForces(const Vector& rExternalForces) { mExternalForces = rExternalForces; }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rResult = {1};
    };

private:
    Vector mInternalForces;
    Vector mExternalForces;
};

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSAtStartOfStageIsEqualToCurrentInternalForcesAtStartMinusInternalForcesAtStart,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadSteppingScheme>();
    element->SetId(1);

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 0.0;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    Vector internal_forces_at_start_of_stage(4);
    internal_forces_at_start_of_stage <<= 1.0, 2.0, 3.0, 4.0; // Change later to new creation function
    Vector external_forces_at_start_of_stage = internal_forces_at_start_of_stage * 2;
    element->SetInternalForces(internal_forces_at_start_of_stage);
    element->SetExternalForces(external_forces_at_start_of_stage);
    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;

    Vector difference{ScalarVector(4, 0.1)};
    Vector current_internal_forces = internal_forces_at_start_of_stage + difference;
    element->SetInternalForces(current_internal_forces);
    scheme.CalculateRHSContribution(*element, actual_right_hand_side, EquationId, CurrentProcessInfo);

    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, difference, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSAtEndOfStageIsEqualToCurrentInternalForcesPlusExternalForcesAtStart,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadSteppingScheme>();
    element->SetId(1);

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 1.0;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    Vector internal_forces_at_start_of_stage(4);
    internal_forces_at_start_of_stage <<= 1.0, 2.0, 3.0, 4.0; // Change later to new creation function

    Vector external_forces_at_start_of_stage = internal_forces_at_start_of_stage * 2;
    element->SetInternalForces(internal_forces_at_start_of_stage);
    element->SetExternalForces(external_forces_at_start_of_stage);
    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;

    Vector difference{ScalarVector(4, 0.1)};
    Vector current_internal_forces = internal_forces_at_start_of_stage + difference;
    element->SetInternalForces(current_internal_forces);
    scheme.CalculateRHSContribution(*element, actual_right_hand_side, EquationId, CurrentProcessInfo);

    auto expected = current_internal_forces + external_forces_at_start_of_stage;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSAtHalfWayIsWeightedSum, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadSteppingScheme>();
    element->SetId(1);

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 0.5;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    Vector internal_forces_at_start_of_stage(4);
    internal_forces_at_start_of_stage <<= -1.0, -2.0, -3.0, -4.0; // Change later to new creation function

    Vector external_forces_at_start_of_stage(4);
    external_forces_at_start_of_stage <<= 5.0, 6.0, 7.0, 8.0;
    element->SetInternalForces(internal_forces_at_start_of_stage);
    element->SetExternalForces(external_forces_at_start_of_stage);
    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;

    Vector current_internal_forces(4);
    current_internal_forces <<= -2.0, -3.0, -4.0, -5.0;
    element->SetInternalForces(current_internal_forces);
    scheme.CalculateRHSContribution(*element, actual_right_hand_side, EquationId, CurrentProcessInfo);

    auto expected = Vector{ScalarVector(4, 1.0)};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected, 1e-6);

    std::vector<std::size_t> expected_equation_ids = {1};
    KRATOS_EXPECT_VECTOR_EQ(EquationId, expected_equation_ids);
}

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeRHSViaLocalSystemAtHalfWayIsWeightedSum, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadSteppingScheme>();
    element->SetId(1);

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[TIME]       = 0.5;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;

    Vector internal_forces_at_start_of_stage(4);
    internal_forces_at_start_of_stage <<= -1.0, -2.0, -3.0, -4.0; // Change later to new creation function

    Vector external_forces_at_start_of_stage(4);
    external_forces_at_start_of_stage <<= 5.0, 6.0, 7.0, 8.0;
    element->SetInternalForces(internal_forces_at_start_of_stage);
    element->SetExternalForces(external_forces_at_start_of_stage);
    scheme.InitializeSolutionStep(model_part, A, Dx, b);
    Vector actual_right_hand_side;

    Vector current_internal_forces(4);
    current_internal_forces <<= -2.0, -3.0, -4.0, -5.0;
    element->SetInternalForces(current_internal_forces);
    Matrix LHS;
    scheme.CalculateSystemContributions(*element, LHS, actual_right_hand_side, EquationId, CurrentProcessInfo);

    auto expected = Vector{ScalarVector(4, 1.0)};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected, 1e-6);

    std::vector<std::size_t> expected_equation_ids = {1};
    KRATOS_EXPECT_VECTOR_EQ(EquationId, expected_equation_ids);
}

class MockConditionForLoadStepping : public Condition
{
public:
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        rRightHandSideVector = Vector{ScalarVector(4, 10.0)};
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rResult = {1};
    };
};

KRATOS_TEST_CASE_IN_SUITE(LoadSteppingSchemeConditionRHSIsScaledWithTime, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    MockConditionForLoadStepping                                           condition;

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    std::vector<std::size_t> EquationId;

    Vector actual_right_hand_side;
    scheme.CalculateRHSContribution(condition, actual_right_hand_side, EquationId, CurrentProcessInfo);


    // This probably should be a parametrized test later
    CurrentProcessInfo[TIME]       = 0.0;
    auto expected = Vector{ZeroVector(4)};
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected, 1e-6);

    CurrentProcessInfo[TIME]       = 0.6;
    scheme.CalculateRHSContribution(condition, actual_right_hand_side, EquationId, CurrentProcessInfo);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, Vector{ScalarVector(4, 6.0)}, 1e-6);


    CurrentProcessInfo[TIME]       = 1.0;
    scheme.CalculateRHSContribution(condition, actual_right_hand_side, EquationId, CurrentProcessInfo);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, Vector{ScalarVector(4, 10.0)}, 1e-6);


    std::vector<std::size_t> expected_equation_ids = {1};
    KRATOS_EXPECT_VECTOR_EQ(EquationId, expected_equation_ids);
}

} // namespace Kratos::Testing