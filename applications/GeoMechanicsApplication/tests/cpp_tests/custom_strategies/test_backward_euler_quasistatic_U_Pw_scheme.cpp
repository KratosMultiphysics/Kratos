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
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;

class BackwardEulerQuasiStaticUPwSchemeTester
{
public:
    Model mModel;

    BackwardEulerQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType> mScheme;

    BackwardEulerQuasiStaticUPwSchemeTester() { CreateValidModelPart(); }

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

        p_node->FastGetSolutionStepValue(DISPLACEMENT, 1) = Kratos::array_1d<double, 3>{7.0, 8.0, 9.0};
        p_node->FastGetSolutionStepValue(VELOCITY, 1) = Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
        p_node->FastGetSolutionStepValue(ACCELERATION, 1) = Kratos::array_1d<double, 3>{4.0, 5.0, 6.0};

        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 1) = 1.0;
        p_node->FastGetSolutionStepValue(WATER_PRESSURE, 0) = 2.0;
    }

    ModelPart& GetModelPart() { return mModel.GetModelPart("dummy"); }
};

KRATOS_TEST_CASE_IN_SUITE(CheckBackwardEulerUPwScheme_ReturnsZeroForValidModelPart, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    BackwardEulerQuasiStaticUPwSchemeTester tester;
    KRATOS_EXPECT_EQ(tester.mScheme.Check(tester.GetModelPart()), 0);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeBackwardEulerUPwScheme_SetsTimeFactors, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    BackwardEulerQuasiStaticUPwSchemeTester tester;

    tester.mScheme.Initialize(tester.GetModelPart());

    // These are the expected numbers according to the SetTimeFactors function
    constexpr double expected_dt_pressure_coefficient = 1.0 / 4.0;
    constexpr double expected_velocity_coefficient    = 1.0 / 4.0;
    KRATOS_EXPECT_TRUE(tester.mScheme.SchemeIsInitialized())

    CompressedMatrix A;
    Vector           Dx;
    Vector           b;
    tester.mScheme.InitializeSolutionStep(tester.GetModelPart(), A, Dx, b); // This is needed to set the time factors

    KRATOS_EXPECT_DOUBLE_EQ(tester.GetModelPart().GetProcessInfo()[DT_PRESSURE_COEFFICIENT],
                            expected_dt_pressure_coefficient);
    KRATOS_EXPECT_DOUBLE_EQ(tester.GetModelPart().GetProcessInfo()[VELOCITY_COEFFICIENT],
                            expected_velocity_coefficient);
}

KRATOS_TEST_CASE_IN_SUITE(BackwardEulerUPwSchemePredict_UpdatesVariablesDerivatives, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    BackwardEulerQuasiStaticUPwSchemeTester tester;

    ModelPart::DofsArrayType dof_set;
    CompressedMatrix         A;
    Vector                   Dx;
    Vector                   b;

    tester.mScheme.InitializeSolutionStep(tester.GetModelPart(), A, Dx, b); // This is needed to set the time factors

    tester.mScheme.Predict(tester.GetModelPart(), dof_set, A, Dx, b);

    constexpr auto expected_dt_water_pressure = 0.25;
    KRATOS_EXPECT_DOUBLE_EQ(tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE, 0),
                            expected_dt_water_pressure);
}

KRATOS_TEST_CASE_IN_SUITE(BackwardEulerUPwSchemeUpdate_DoesNotUpdateFixedScalarVariable,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    BackwardEulerQuasiStaticUPwSchemeTester tester;

    ModelPart::DofsArrayType dof_set;
    CompressedMatrix         A;
    Vector                   Dx;
    Vector                   b;

    tester.mScheme.InitializeSolutionStep(tester.GetModelPart(), A, Dx, b); // This is needed to set the time factors

    tester.GetModelPart().Nodes()[0].Fix(DT_WATER_PRESSURE);

    tester.mScheme.Update(tester.GetModelPart(), dof_set, A, Dx, b);

    // should be the same as the original
    const auto expected_dt_water_pressure = 0.0;

    const auto actual_dt_water_pressure =
        tester.GetModelPart().Nodes()[0].FastGetSolutionStepValue(DT_WATER_PRESSURE, 0);
    KRATOS_EXPECT_DOUBLE_EQ(actual_dt_water_pressure, expected_dt_water_pressure);
}

} // namespace Kratos::Testing
