// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include "custom_processes/calculate_total_motion_process.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateTotalMotionProcessDisplacement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);
    r_model_part.AddNodalSolutionStepVariable(INCREMENTAL_DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);

    auto p_node_1 = r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    std::vector<intrusive_ptr<Node>> nodes = {p_node_1, p_node_2};

    p_node_1->FastGetSolutionStepValue(TOTAL_DISPLACEMENT, 0) = Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node_1->FastGetSolutionStepValue(INCREMENTAL_DISPLACEMENT, 0) =
        Kratos::array_1d<double, 3>{-4.0, -5.0, -6.0};

    p_node_2->FastGetSolutionStepValue(TOTAL_DISPLACEMENT, 0) = Kratos::array_1d<double, 3>{7.0, 12.0, 9.0};
    p_node_2->FastGetSolutionStepValue(INCREMENTAL_DISPLACEMENT, 0) =
        Kratos::array_1d<double, 3>{6.0, 11.0, 8.0};

    auto process =
        CalculateTotalMotionProcess(r_model_part, Parameters(R"({"variable_name": "DISPLACEMENT"})"));

    process.Execute();

    const auto expected_total_displacement_1 = Kratos::array_1d<double, 3>{-3.0, -3.0, -3.0};
    const auto expected_total_displacement_2 = Kratos::array_1d<double, 3>{13.0, 23.0, 17.0};

    KRATOS_CHECK_VECTOR_NEAR(p_node_1->FastGetSolutionStepValue(TOTAL_DISPLACEMENT),
                             expected_total_displacement_1, Defaults::absolute_tolerance)
    KRATOS_CHECK_VECTOR_NEAR(p_node_2->FastGetSolutionStepValue(TOTAL_DISPLACEMENT),
                             expected_total_displacement_2, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateTotalMotionProcessRotation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);
    r_model_part.AddNodalSolutionStepVariable(INCREMENTAL_ROTATION);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_ROTATION);

    auto p_node_1 = r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    std::vector<intrusive_ptr<Node>> nodes = {p_node_1, p_node_2};

    p_node_1->FastGetSolutionStepValue(TOTAL_ROTATION, 0) = Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node_1->FastGetSolutionStepValue(INCREMENTAL_ROTATION, 0) =
        Kratos::array_1d<double, 3>{-4.0, -5.0, -6.0};

    p_node_2->FastGetSolutionStepValue(TOTAL_ROTATION, 0) = Kratos::array_1d<double, 3>{7.0, 12.0, 9.0};
    p_node_2->FastGetSolutionStepValue(INCREMENTAL_ROTATION, 0) = Kratos::array_1d<double, 3>{6.0, 11.0, 8.0};

    auto process = CalculateTotalMotionProcess(r_model_part, Parameters(R"({"variable_name": "ROTATION"})"));

    process.Execute();

    const auto expected_total_rotation_1 = Kratos::array_1d<double, 3>{-3.0, -3.0, -3.0};
    const auto expected_total_rotation_2 = Kratos::array_1d<double, 3>{13.0, 23.0, 17.0};

    KRATOS_CHECK_VECTOR_NEAR(p_node_1->FastGetSolutionStepValue(TOTAL_ROTATION),
                             expected_total_rotation_1, Defaults::absolute_tolerance)
    KRATOS_CHECK_VECTOR_NEAR(p_node_2->FastGetSolutionStepValue(TOTAL_ROTATION),
                             expected_total_rotation_2, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateTotalMotionProcessUndefined, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CalculateTotalMotionProcess(r_model_part, Parameters(R"({"variable_name": "NOTHING"})")),
        "Invalid variable name: NOTHING. Expected DISPLACEMENT or ROTATION.");
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoCalculateTotalMotionProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model                             model;
    auto&                             r_model_part = model.CreateModelPart("dummy", 2);
    const CalculateTotalMotionProcess process(r_model_part, {R"({"variable_name": "ROTATION"})"});

    KRATOS_EXPECT_EQ(process.Info(), "CalculateTotalMotionProcess");
}
} // namespace Kratos::Testing
