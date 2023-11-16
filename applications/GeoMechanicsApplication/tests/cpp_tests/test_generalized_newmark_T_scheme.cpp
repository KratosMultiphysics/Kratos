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

GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> CreateValidScheme()
{
    return GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType>(0.75);
}

ModelPart& CreateValidModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("dummy", 2);
    result.AddNodalSolutionStepVariable(TEMPERATURE);
    result.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = result.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    return result;
}


KRATOS_TEST_CASE_IN_SUITE(CheckBackwardEulerQuasistaticTScheme_WithAllNecessaryParts_Returns0,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto scheme = CreateValidScheme();
    const auto& model_part = CreateValidModelPart(model);

    KRATOS_EXPECT_EQ(scheme.Check(model_part), 0);
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidTheta_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    constexpr double invalid_theta = -2.0;
    GeneralizedNewmarkTScheme<SparseSpaceType, LocalSpaceType> scheme(invalid_theta);

    Model model;
    const auto& model_part = CreateValidModelPart(model);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "Theta has an invalid value")
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidBufferSize_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();

    Model model;
    constexpr auto invalid_buffer_size = ModelPart::IndexType{1};
    auto& model_part = CreateValidModelPart(model);
    model_part.SetBufferSize(invalid_buffer_size);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "insufficient buffer size. Buffer size should be greater than or equal to "
        "2. Current size is 1")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "missing TEMPERATURE dof on node ")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingDtTemperatureSolutionStepVariable_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "DT_TEMPERATURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingTemperatureSolutionStepVariable_CheckBackwardEulerQuasistaticTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "TEMPERATURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ThermalSchemeUpdate_SetsDtTemperature, KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);

    constexpr double current_temperature = 10.0;
    constexpr double previous_temperature = 5.0;
    constexpr double previous_dt_temperature = 3.0;
    constexpr double delta_time = 2.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    Node& node = model_part.Nodes()[0];
    node.FastGetSolutionStepValue(TEMPERATURE) = current_temperature;
    node.FastGetSolutionStepValue(TEMPERATURE, 1) = previous_temperature;
    node.FastGetSolutionStepValue(DT_TEMPERATURE, 1) = previous_dt_temperature;

    // This is the expected value as calculated by the UpdateVariablesDerivatives
    constexpr double expected_value = 7.0/3.0;

    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;

    scheme.Initialize(model_part);
    scheme.Predict(model_part, dof_set, A, Dx, b);

    KRATOS_EXPECT_DOUBLE_EQ(node.FastGetSolutionStepValue(DT_TEMPERATURE), expected_value);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    auto scheme = CreateValidScheme();
    Model model;
    ModelPart& model_part = CreateValidModelPart(model);
    constexpr double delta_time = 3.0;
    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;

    scheme.Initialize(model_part);

    KRATOS_EXPECT_TRUE(scheme.SchemeIsInitialized())
    constexpr double theta = 0.75;
    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[DT_TEMPERATURE_COEFFICIENT],
                            1.0 / (theta * delta_time));
}

} // namespace Kratos::Test