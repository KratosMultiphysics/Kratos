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

#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(BackwardEulerPwScheme_UpdatesVariablesDerivatives_WhenPredictIsCalled,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme;
    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);

    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);

    constexpr double current_pressure = 1.0;
    constexpr double previous_pressure = 0.0;
    constexpr double delta_time = 4.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->FastGetSolutionStepValue(WATER_PRESSURE, 0) = current_pressure;
    p_node->FastGetSolutionStepValue(WATER_PRESSURE, 1) = previous_pressure;

    KRATOS_EXPECT_DOUBLE_EQ(p_node->FastGetSolutionStepValue(DT_WATER_PRESSURE, 0), 0.0);

    scheme.Initialize(model_part);
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;
    scheme.Predict(model_part, dof_set, A, Dx, b);

    constexpr double expected_dt_temperature = 0.25;
    KRATOS_EXPECT_DOUBLE_EQ(p_node->FastGetSolutionStepValue(DT_WATER_PRESSURE, 0),
                            expected_dt_temperature);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeBackwardEulerPwScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme;
    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);

    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);

    constexpr double delta_time = 3.0;
    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;

    scheme.Initialize(model_part);

    KRATOS_EXPECT_TRUE(scheme.SchemeIsInitialized())
    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[DT_PRESSURE_COEFFICIENT],
                            1.0 / delta_time);
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckBackwardEulerPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "missing WATER_PRESSURE dof on node ")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingDtWaterPressureSolutionStepVariable_CheckBackwardEulerPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(WATER_PRESSURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "DT_WATER_PRESSURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingWaterPressureSolutionStepVariable_CheckBackwardEulerPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "WATER_PRESSURE variable is not allocated for node 0")
}

} // namespace Kratos::Testing
