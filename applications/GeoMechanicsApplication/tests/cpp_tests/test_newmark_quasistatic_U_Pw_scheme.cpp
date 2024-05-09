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

#include "testing/testing.h"

#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "spaces/ublas_space.h"
#include "test_utilities/spy_element.h"
#include "test_utilities/spy_condition.h"

namespace Kratos::Testing {

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

class NewmarkQuasistaticUPwSchemeTester {
public:
    Model mModel;
    NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType> mScheme =
        CreateValidScheme();

    NewmarkQuasistaticUPwSchemeTester()
    {
        CreateValidModelPart();
    }

    NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType> CreateValidScheme() const
    {
        return NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType>(0.25, 0.5, 0.75);
    }

    void CreateValidModelPart()
    {
        auto& result = mModel.CreateModelPart("dummy", 2);
        result.AddNodalSolutionStepVariable(VELOCITY);
        result.AddNodalSolutionStepVariable(ACCELERATION);
        result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
        result.AddNodalSolutionStepVariable(DISPLACEMENT);
        result.AddNodalSolutionStepVariable(WATER_PRESSURE);

        auto p_node = result.CreateNewNode(0, 0.0, 0.0, 0.0);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->AddDof(WATER_PRESSURE);
        result.GetProcessInfo()[DELTA_TIME] = 4.0;

        p_node->FastGetSolutionStepValue(DISPLACEMENT, 1) =
            Kratos::array_1d<double, 3>{7.0, 8.0, 9.0};
        p_node->FastGetSolutionStepValue(VELOCITY, 1) =
            Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 1) =
            Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};

        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 1) = 1.0;
        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 0) = 2.0;
    }

    ModelPart& GetModelPart()
    {
        return mModel.GetModelPart("dummy");
    }
};

KRATOS_TEST_CASE_IN_SUITE(CheckNewmarkUPwScheme_ReturnsZeroForValidScheme, KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;
    KRATOS_EXPECT_EQ(tester.mScheme.Check(tester.GetModelPart()), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeNewmarkUPwScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart());

    // These are the expected numbers according to the SetTimeFactors function
    constexpr double expected_dt_pressure_coefficient = 1.0 / 3.0;
    constexpr double expected_velocity_coefficient = 0.5;
    KRATOS_EXPECT_TRUE(tester.mScheme.SchemeIsInitialized())
    KRATOS_EXPECT_DOUBLE_EQ(tester.GetModelPart().GetProcessInfo()[DT_PRESSURE_COEFFICIENT],
                            expected_dt_pressure_coefficient);
    KRATOS_EXPECT_DOUBLE_EQ(tester.GetModelPart().GetProcessInfo()[VELOCITY_COEFFICIENT],
                            expected_velocity_coefficient);
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkUPwSchemePredict_UpdatesVariablesDerivatives, KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // These expected numbers result from the calculations in UpdateVariablesDerivatives
    const auto expected_acceleration = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
    const auto expected_velocity = Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
    constexpr auto expected_dt_water_pressure = 1.0 / 3.0;

    const auto actual_acceleration = tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(ACCELERATION, 0);
    KRATOS_EXPECT_EQ(actual_acceleration.size(), expected_acceleration.size());
    for (std::size_t i = 0; i < actual_acceleration.size(); ++i)
        KRATOS_EXPECT_DOUBLE_EQ(actual_acceleration[i], expected_acceleration[i]);

    const auto actual_velocity = tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(VELOCITY, 0);
    for (std::size_t i = 0; i < actual_velocity.size(); ++i)
        KRATOS_EXPECT_DOUBLE_EQ(actual_velocity[i], expected_velocity[i]);

    KRATOS_EXPECT_DOUBLE_EQ(
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE, 0),
        expected_dt_water_pressure);
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckNewmarkUPwScheme_Throws, KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticUPwSchemeTester tester;
    auto& model_part = tester.mModel.CreateModelPart("MissingAcceleration", 2);
    model_part.AddNodalSolutionStepVariable(VELOCITY); // Missing ACCELERATION
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(DISPLACEMENT_X);
    p_node->AddDof(DISPLACEMENT_Y);
    p_node->AddDof(DISPLACEMENT_Z);
    p_node->AddDof(WATER_PRESSURE);
    model_part.GetProcessInfo()[DELTA_TIME] = 4.0;

    p_node->FastGetSolutionStepValue(VELOCITY, 1) =
        Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};
    p_node->FastGetSolutionStepValue(DISPLACEMENT, 1) =
        Kratos::array_1d<double, 3>{0.0, 0.0, 0.0};

    p_node->FastGetSolutionStepValue(WATER_PRESSURE, 1) = 1.0;
    p_node->FastGetSolutionStepValue(WATER_PRESSURE, 0) = 2.0;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(tester.mScheme.Check(model_part), "ACCELERATION variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidBeta_CheckNewmarkUPwScheme_Throws, KratosGeoMechanicsFastSuite)
{
    constexpr double invalid_beta = -2.3;
    using SchemeType = NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType>;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        SchemeType scheme(invalid_beta, 0.5, 0.75),
        "Beta must be larger than zero, but got -2.3")
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidGamma_CheckNewmarkUPwScheme_Throws, KratosGeoMechanicsFastSuite)
{
    constexpr double invalid_gamma = -2.5;
    using SchemeType = NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType>;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        SchemeType scheme(0.25, invalid_gamma, 0.75),
        "Gamma must be larger than zero, but got -2.5")
}



} // namespace Kratos::Testing