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

#include "custom_strategies/schemes/generalized_newmark_T_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing {

ModelPart& CreateValidTemperatureModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy", 2);
    result.AddNodalSolutionStepVariable(TEMPERATURE);
    result.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = result.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    return result;
}

KRATOS_TEST_CASE_IN_SUITE(CheckNewmarkTScheme_WithAllNecessaryParts_Returns0, KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    const auto& model_part = CreateValidTemperatureModelPart(model);

    KRATOS_EXPECT_EQ(scheme.Check(model_part), 0);
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckNewmarkTScheme_Throws, KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "missing TEMPERATURE dof on node ")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingDtTemperatureSolutionStepVariable_CheckNewmarkTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "DT_TEMPERATURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingTemperatureSolutionStepVariable_CheckNewmarkTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "TEMPERATURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(NewmarkTSchemeUpdate_SetsDtTemperature, KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    ModelPart& model_part = CreateValidTemperatureModelPart(model);

    constexpr double current_temperature = 10.0;
    constexpr double previous_temperature = 5.0;
    constexpr double previous_dt_temperature = 3.0;
    constexpr double delta_time = 2.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    Node& node = model_part.Nodes()[0];
    node.FastGetSolutionStepValue(TEMPERATURE, 0) = current_temperature;
    node.FastGetSolutionStepValue(TEMPERATURE, 1) = previous_temperature;
    node.FastGetSolutionStepValue(DT_TEMPERATURE, 1) = previous_dt_temperature;

    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    scheme.Initialize(model_part);
    scheme.Predict(model_part, dof_set, A, Dx, b);

    // This is the expected value as calculated by the UpdateVariablesDerivatives
    KRATOS_EXPECT_DOUBLE_EQ(node.FastGetSolutionStepValue(DT_TEMPERATURE, 0), 7.0 / 3.0);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeNewmarkTScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    constexpr double theta = 0.75;
    GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> scheme(theta);

    Model model;
    ModelPart& model_part = CreateValidTemperatureModelPart(model);
    constexpr double delta_time = 3.0;
    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;

    scheme.Initialize(model_part);

    KRATOS_EXPECT_TRUE(scheme.SchemeIsInitialized())
    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[DT_TEMPERATURE_COEFFICIENT],
                            1.0 / (theta * delta_time));
}

} // namespace Kratos::Testing