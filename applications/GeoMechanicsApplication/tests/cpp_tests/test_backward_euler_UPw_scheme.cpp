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

#include "containers/model.h"
#include "custom_strategies/schemes/backward_euler_quasistatic_U_Pw_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

class BackwardEulerUPwSchemeTester
{
public:
    Model mModel;
    BackwardEulerQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType> mScheme;

    BackwardEulerUPwSchemeTester() { CreateValidModelPart(); }

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
            Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 1) =
            Kratos::array_1d<double, 3>{4.0, 5.0, 6.0};

        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 1) = 1.0;
        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 0) = 2.0;
    }

    ModelPart& GetModelPart() { return mModel.GetModelPart("dummy"); }
};

KRATOS_TEST_CASE_IN_SUITE(CheckBackwardEulerUPwScheme_ReturnsZeroForValidModelPart,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerUPwSchemeTester tester;
    KRATOS_EXPECT_EQ(tester.mScheme.Check(tester.GetModelPart()), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeBackwardEulerUPwScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    BackwardEulerUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart());

    // These are the expected numbers according to the SetTimeFactors function
    constexpr double expected_dt_pressure_coefficient = 1.0 / 4.0;
    constexpr double expected_velocity_coefficient = 1.0 / 4.0;
    KRATOS_EXPECT_TRUE(tester.mScheme.SchemeIsInitialized())
    KRATOS_EXPECT_DOUBLE_EQ(tester.GetModelPart().GetProcessInfo()[DT_PRESSURE_COEFFICIENT],
                            expected_dt_pressure_coefficient);
    KRATOS_EXPECT_DOUBLE_EQ(tester.GetModelPart().GetProcessInfo()[VELOCITY_COEFFICIENT],
                            expected_velocity_coefficient);
}

KRATOS_TEST_CASE_IN_SUITE(BackwardEulerUPwSchemePredict_UpdatesVariablesDerivatives,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart()); // This is needed to set the time factors
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    // These expected numbers result from the calculations in UpdateVariablesDerivatives
    const auto expected_acceleration = Kratos::array_1d<double, 3>{-0.6875, -1.0, -1.3125};
    const auto expected_velocity = Kratos::array_1d<double, 3>{-1.75, -2.0, -2.25};
    constexpr auto expected_dt_water_pressure = 0.25;

    const auto actual_acceleration =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(ACCELERATION, 0);
    const auto actual_velocity =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(VELOCITY, 0);

    constexpr auto absolute_tolerance = 1.0e-6;
    KRATOS_EXPECT_VECTOR_NEAR(expected_acceleration, actual_acceleration, absolute_tolerance)
    KRATOS_EXPECT_VECTOR_NEAR(expected_velocity, actual_velocity, absolute_tolerance)

    KRATOS_EXPECT_DOUBLE_EQ(
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE, 0),
        expected_dt_water_pressure);
}

} // namespace Kratos::Testing
