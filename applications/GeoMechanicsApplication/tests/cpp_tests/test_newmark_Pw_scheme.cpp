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

#include "custom_strategies/schemes/newmark_quasistatic_Pw_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing {

namespace {
NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> CreateValidScheme()
{
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> result(0.75);
    return result;
}

ModelPart& CreateValidModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy", 2);
    result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    auto p_node = result.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(WATER_PRESSURE);

    return result;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(CheckNewmarkPwScheme_WithAllNecessaryParts_Returns0,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto scheme = CreateValidScheme();
    const auto& model_part = CreateValidModelPart(model);

    KRATOS_EXPECT_EQ(scheme.Check(model_part), 0);
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidTheta_CheckNewmarkPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    constexpr int invalid_theta = -2;
    using SchemeType = NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(SchemeType scheme(invalid_theta),
                                      "Theta must be larger than zero, but got -2")
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidBufferSize_CheckNewmarkPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    constexpr int invalid_buffer_size = 1;
    auto& model_part = model.CreateModelPart("dummy", invalid_buffer_size);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(WATER_PRESSURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "insufficient buffer size. Buffer size should be greater than or equal to "
        "2. Current size is 1")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckNewmarkPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "missing WATER_PRESSURE dof on node ")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingDtWaterPressureSolutionStepVariable_CheckNewmarkPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(WATER_PRESSURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "DT_WATER_PRESSURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingWaterPressureSolutionStepVariable_CheckNewmarkPwScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(0.75);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "WATER_PRESSURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkPwSchemeUpdate_SetsDtPressure, KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);

    constexpr double current_pressure = 10.0;
    constexpr double previous_temperature = 5.0;
    constexpr double previous_dt_temperature = 3.0;
    constexpr double delta_time = 2.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    Node& node = model_part.Nodes()[0];
    node.FastGetSolutionStepValue(WATER_PRESSURE) = current_pressure;
    node.FastGetSolutionStepValue(WATER_PRESSURE, 1) = previous_temperature;
    node.FastGetSolutionStepValue(DT_WATER_PRESSURE, 1) = previous_dt_temperature;

    // This is the expected value as calculated by the UpdateVariablesDerivatives
    constexpr double expected_value = 7.0/3.0;

    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    scheme.Initialize(model_part);
    scheme.Predict(model_part, dof_set, A, Dx, b);

    KRATOS_EXPECT_DOUBLE_EQ(node.FastGetSolutionStepValue(DT_WATER_PRESSURE), expected_value);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeNewmarkPwScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);
    constexpr double delta_time = 3.0;
    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;

    scheme.Initialize(model_part);

    KRATOS_EXPECT_TRUE(scheme.SchemeIsInitialized())
    constexpr double theta = 0.75;
    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[DT_PRESSURE_COEFFICIENT],
                            1.0 / (theta * delta_time));
}

} // namespace Kratos::Testing