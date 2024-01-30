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

KRATOS_TEST_CASE_IN_SUITE(CheckNewmarkPwScheme_WithAllNecessaryParts_Returns0, KratosGeoMechanicsFastSuite)
{
    Model model;
    constexpr double theta = 0.75;
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(theta);
    const auto& model_part = CreateValidModelPart(model);

    KRATOS_EXPECT_EQ(scheme.Check(model_part), 0);
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckNewmarkPwScheme_Throws, KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

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
    constexpr double theta = 0.75;
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

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
    constexpr double theta = 0.75;
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "WATER_PRESSURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkPwSchemeUpdate_SetsDtPressure, KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    ModelPart& model_part = CreateValidModelPart(model);

    constexpr double current_pressure = 10.0;
    constexpr double previous_pressure = 5.0;
    constexpr double previous_dt_pressure = 3.0;
    constexpr double delta_time = 2.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    Node& node = model_part.Nodes()[0];
    node.FastGetSolutionStepValue(WATER_PRESSURE, 0) = current_pressure;
    node.FastGetSolutionStepValue(WATER_PRESSURE, 1) = previous_pressure;
    node.FastGetSolutionStepValue(DT_WATER_PRESSURE, 1) = previous_dt_pressure;

    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    scheme.Initialize(model_part);
    scheme.Predict(model_part, dof_set, A, Dx, b);

    // This is the expected value as calculated by the UpdateVariablesDerivatives
    KRATOS_EXPECT_DOUBLE_EQ(node.FastGetSolutionStepValue(DT_WATER_PRESSURE, 0), 7.0 / 3.0);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeNewmarkPwScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    ModelPart& model_part = CreateValidModelPart(model);
    constexpr double delta_time = 3.0;
    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;

    scheme.Initialize(model_part);

    KRATOS_EXPECT_TRUE(scheme.SchemeIsInitialized())
    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[DT_PRESSURE_COEFFICIENT],
                            1.0 / (theta * delta_time));
}

} // namespace Kratos::Testing