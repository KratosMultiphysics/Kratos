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

#include "custom_strategies/schemes/backward_euler_T_scheme.hpp"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

namespace Kratos::Testing {

using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

KRATOS_TEST_CASE_IN_SUITE(BackwardEulerTScheme_UpdatesVariablesDerivatives_WhenPredictIsCalled,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerTScheme<SparseSpaceType, LocalSpaceType> scheme;
    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);

    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);

    constexpr double current_temperature = 1.0;
    constexpr double previous_temperature = 0.0;
    constexpr double delta_time = 4.0;

    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->FastGetSolutionStepValue(TEMPERATURE) = current_temperature;
    p_node->FastGetSolutionStepValue(TEMPERATURE, 1) = previous_temperature;

    KRATOS_EXPECT_DOUBLE_EQ(p_node->FastGetSolutionStepValue(DT_TEMPERATURE), 0.0);

    scheme.Initialize(model_part);
    ModelPart::DofsArrayType dof_set;
    CompressedMatrix A;
    Vector Dx;
    Vector b;
    scheme.Predict(model_part, dof_set, A, Dx, b);

    constexpr double expected_dt_temperature = 0.25;
    KRATOS_EXPECT_DOUBLE_EQ(p_node->FastGetSolutionStepValue(DT_TEMPERATURE),
                            expected_dt_temperature);
}

KRATOS_TEST_CASE_IN_SUITE(InitializeBackwardEulerTScheme_SetsTimeFactors, KratosGeoMechanicsFastSuite)
{
    BackwardEulerTScheme<SparseSpaceType, LocalSpaceType> scheme;
    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);

    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);

    constexpr double delta_time = 3.0;
    model_part.GetProcessInfo()[DELTA_TIME] = delta_time;

    scheme.Initialize(model_part);

    KRATOS_EXPECT_TRUE(scheme.SchemeIsInitialized())
    KRATOS_EXPECT_DOUBLE_EQ(model_part.GetProcessInfo()[DT_TEMPERATURE_COEFFICIENT],
                            1.0 / (delta_time));
}

KRATOS_TEST_CASE_IN_SUITE(ForInvalidBufferSize_CheckBackwardEulerTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerTScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    constexpr int invalid_buffer_size = 1;
    auto& model_part = model.CreateModelPart("dummy", invalid_buffer_size);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "insufficient buffer size. Buffer size should be greater than or equal to "
        "2. Current size is 1")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingNodalDof_CheckBackwardEulerTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerTScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(scheme.Check(model_part),
                                      "missing TEMPERATURE dof on node ")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingDtTemperatureSolutionStepVariable_CheckBackwardEulerTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerTScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    p_node->AddDof(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "DT_TEMPERATURE variable is not allocated for node 0")
}

KRATOS_TEST_CASE_IN_SUITE(ForMissingTemperatureSolutionStepVariable_CheckBackwardEulerTScheme_Throws,
                          KratosGeoMechanicsFastSuite)
{
    BackwardEulerTScheme<SparseSpaceType, LocalSpaceType> scheme;

    Model model;
    auto& model_part = model.CreateModelPart("dummy", 2);
    model_part.AddNodalSolutionStepVariable(DT_TEMPERATURE);
    auto p_node = model_part.CreateNewNode(0, 0.0, 0.0, 0.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        scheme.Check(model_part),
        "TEMPERATURE variable is not allocated for node 0")
}

} // namespace Kratos::Testing