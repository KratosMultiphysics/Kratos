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
#include "custom_utilities/ublas_utilities.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace
{

using namespace Kratos;

class MockElementForLoadStepping final : public Element
{
public:
    void Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo&) override
    {
        if (rVariable == INTERNAL_FORCES_VECTOR) Output = mInternalForces;
        if (rVariable == EXTERNAL_FORCES_VECTOR) Output = mExternalForces;
    }

    using Element::Calculate;

    void SetInternalForces(const Vector& rInternalForces) { mInternalForces = rInternalForces; }

    void SetExternalForces(const Vector& rExternalForces) { mExternalForces = rExternalForces; }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override
    {
        rResult = {3, 4};
    }

private:
    Vector mInternalForces;
    Vector mExternalForces;
};

class MockConditionForLoadStepping final : public Condition
{
public:
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        rRightHandSideVector = Vector(4, 10.0);
    }

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override
    {
        rResult = {1, 2};
    }
};

} // namespace

namespace Kratos::Testing
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

class LoadSteppingSchemeElementRightHandSideScaling
    : public ::testing::TestWithParam<std::tuple<double, Vector>>
{
};

TEST_P(LoadSteppingSchemeElementRightHandSideScaling, RightHandSideIsCalculatedBasedOnTime)
{
    // Arrange
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    auto element = Kratos::make_intrusive<MockElementForLoadStepping>();
    element->SetInternalForces(UblasUtilities::CreateVector({1.0, 2.0, 3.0, 4.0}));
    element->SetExternalForces(UblasUtilities::CreateVector({5.0, 6.0, 7.0, 8.0}));

    ProcessInfo CurrentProcessInfo;
    const auto& [current_time, expected_right_hand_side] = GetParam();
    CurrentProcessInfo[TIME]                             = current_time;
    CurrentProcessInfo[START_TIME]                       = 0.0;
    CurrentProcessInfo[END_TIME]                         = 1.0;
    std::vector<std::size_t> equation_ids;

    Model model;
    auto& model_part = model.CreateModelPart("Main");
    model_part.AddElement(element);
    CompressedMatrix A;
    Vector           Dx;
    Vector           b;
    scheme.InitializeSolutionStep(model_part, A, Dx, b);

    element->SetInternalForces(UblasUtilities::CreateVector({2.0, 3.0, 4.0, 5.0}));

    // Act & Assert
    const std::vector expected_equation_ids = {3, 4};
    Vector            actual_right_hand_side;
    scheme.CalculateRHSContribution(*element, actual_right_hand_side, equation_ids, CurrentProcessInfo);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance);
    KRATOS_EXPECT_VECTOR_EQ(equation_ids, expected_equation_ids);

    Vector actual_right_hand_side_via_system_contribution;
    Matrix dummy_matrix;
    scheme.CalculateSystemContributions(*element, dummy_matrix, actual_right_hand_side_via_system_contribution,
                                        equation_ids, CurrentProcessInfo);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance);
    KRATOS_EXPECT_VECTOR_EQ(equation_ids, expected_equation_ids);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         LoadSteppingSchemeElementRightHandSideScaling,
                         ::testing::Values(std::make_tuple(0.0, Vector(4, -1.0)),
                                           std::make_tuple(0.5, Vector(4, 1.0)),
                                           std::make_tuple(0.8, Vector(4, 2.2)),
                                           std::make_tuple(1.0, Vector(4, 3.0))));

class LoadSteppingSchemeConditionRightHandSideScaling
    : public ::testing::TestWithParam<std::tuple<double, Vector>>
{
};

TEST_P(LoadSteppingSchemeConditionRightHandSideScaling, RightHandSideIsScaledBasedOnTime)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    MockConditionForLoadStepping                        condition;
    const auto& [current_time, expected_right_hand_side] = GetParam();

    ProcessInfo CurrentProcessInfo;
    CurrentProcessInfo[START_TIME] = 0.0;
    CurrentProcessInfo[END_TIME]   = 1.0;
    CurrentProcessInfo[TIME]       = current_time;
    std::vector<std::size_t> actual_equation_ids;

    Vector actual_right_hand_side;
    scheme.CalculateRHSContribution(condition, actual_right_hand_side, actual_equation_ids, CurrentProcessInfo);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side, expected_right_hand_side, Defaults::relative_tolerance);

    Vector actual_right_hand_side_via_system_contribution;
    Matrix dummy_matrix;
    scheme.CalculateSystemContributions(condition, dummy_matrix, actual_right_hand_side_via_system_contribution,
                                        actual_equation_ids, CurrentProcessInfo);
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(actual_right_hand_side_via_system_contribution,
                                       expected_right_hand_side, Defaults::relative_tolerance);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         LoadSteppingSchemeConditionRightHandSideScaling,
                         ::testing::Values(std::make_tuple(0.0, Vector(4, 0.0)),
                                           std::make_tuple(0.6, Vector(4, 6.0)),
                                           std::make_tuple(1.0, Vector(4, 10.0))));

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, LoadSteppingScheme_CalculateRHSContribution_ReturnsCorrectConditionEquationIds)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    MockConditionForLoadStepping                        condition;

    std::vector<std::size_t> actual_equation_ids;
    Vector                   actual_right_hand_side;
    scheme.CalculateRHSContribution(condition, actual_right_hand_side, actual_equation_ids, ProcessInfo{});

    const std::vector expected_equation_ids = {1, 2};
    KRATOS_EXPECT_VECTOR_EQ(actual_equation_ids, expected_equation_ids);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       LoadSteppingScheme_CalculateSystemContributions_ReturnsCorrectConditionEquationIds)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    MockConditionForLoadStepping                        condition;

    std::vector<std::size_t> actual_equation_ids;
    Vector                   actual_right_hand_side;
    Matrix                   dummy_matrix;
    scheme.CalculateSystemContributions(condition, dummy_matrix, actual_right_hand_side,
                                        actual_equation_ids, ProcessInfo{});

    const std::vector expected_equation_ids = {1, 2};
    KRATOS_EXPECT_VECTOR_EQ(actual_equation_ids, expected_equation_ids);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       LoadSteppingScheme_CalculateElementRHSThrows_WhenCallingBeforeInitializeSolutionStep)
{
    LoadSteppingScheme<SparseSpaceType, LocalSpaceType> scheme;
    MockElementForLoadStepping                          element;

    std::vector<std::size_t> actual_equation_ids;
    Vector                   actual_right_hand_side;
    const std::string        expecteded_error_message =
        "The load stepping scheme is not initialized properly. Make sure all initializations are "
        "done before calculating system contributions.";
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.CalculateRHSContribution(element, actual_right_hand_side, actual_equation_ids, ProcessInfo{}),
        expecteded_error_message);
}

} // namespace Kratos::Testing