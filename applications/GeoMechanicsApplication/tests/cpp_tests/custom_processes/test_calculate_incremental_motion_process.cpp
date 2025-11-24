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

#include "custom_processes/calculate_incremental_motion_process.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos::Testing
{

using namespace Kratos;

KRATOS_TEST_CASE_IN_SUITE(CalculateIncrementalMotionProcessDisplacement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // set up model part
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(INCREMENTAL_DISPLACEMENT);

    auto p_node_1 = r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    std::vector<intrusive_ptr<Node>> nodes = {p_node_1, p_node_2};

    for (auto p_node : nodes) {
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
    }

    // set up displacement values
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT, 0) = Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT, 1) = Kratos::array_1d<double, 3>{4.0, 5.0, 6.0};

    p_node_2->FastGetSolutionStepValue(DISPLACEMENT, 0) = Kratos::array_1d<double, 3>{7.0, 12.0, 9.0};
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT, 1) = Kratos::array_1d<double, 3>{11.0, 8.0, 14.0};

    // initialize process to be tested
    auto process = CalculateIncrementalMotionProcess(
        r_model_part, Parameters(R"({"variable_name": "DISPLACEMENT"})"));
    process.Execute();

    const auto expected_incremental_displacement_1 = Kratos::array_1d<double, 3>{-3.0, -3.0, -3.0};
    const auto expected_incremental_displacement_2 = Kratos::array_1d<double, 3>{-4.0, 4.0, -5.0};

    KRATOS_CHECK_VECTOR_NEAR(p_node_1->FastGetSolutionStepValue(INCREMENTAL_DISPLACEMENT),
                             expected_incremental_displacement_1, Defaults::absolute_tolerance)
    KRATOS_CHECK_VECTOR_NEAR(p_node_2->FastGetSolutionStepValue(INCREMENTAL_DISPLACEMENT),
                             expected_incremental_displacement_2, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateIncrementalMotionProcessRotation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // set up model part
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);
    r_model_part.AddNodalSolutionStepVariable(ROTATION);
    r_model_part.AddNodalSolutionStepVariable(INCREMENTAL_ROTATION);

    auto p_node_1 = r_model_part.CreateNewNode(0, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);

    std::vector<intrusive_ptr<Node>> nodes = {p_node_1, p_node_2};

    for (auto p_node : nodes) {
        p_node->AddDof(ROTATION_X);
        p_node->AddDof(ROTATION_Y);
        p_node->AddDof(ROTATION_Z);
    }

    p_node_1->FastGetSolutionStepValue(ROTATION, 0) = Kratos::array_1d<double, 3>{1.0, 2.0, 3.0};
    p_node_1->FastGetSolutionStepValue(ROTATION, 1) = Kratos::array_1d<double, 3>{4.0, 5.0, 6.0};

    p_node_2->FastGetSolutionStepValue(ROTATION, 0) = Kratos::array_1d<double, 3>{7.0, 12.0, 9.0};
    p_node_2->FastGetSolutionStepValue(ROTATION, 1) = Kratos::array_1d<double, 3>{11.0, 8.0, 14.0};

    // initialize process to be tested
    auto process =
        CalculateIncrementalMotionProcess(r_model_part, Parameters(R"({"variable_name": "ROTATION"})"));
    process.Execute();

    const auto expected_incremental_rotation_1 = Kratos::array_1d<double, 3>{-3.0, -3.0, -3.0};
    const auto expected_incremental_rotation_2 = Kratos::array_1d<double, 3>{-4.0, 4.0, -5.0};

    KRATOS_CHECK_VECTOR_NEAR(p_node_1->FastGetSolutionStepValue(INCREMENTAL_ROTATION),
                             expected_incremental_rotation_1, Defaults::absolute_tolerance)
    KRATOS_CHECK_VECTOR_NEAR(p_node_2->FastGetSolutionStepValue(INCREMENTAL_ROTATION),
                             expected_incremental_rotation_2, Defaults::absolute_tolerance)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateIncrementalMotionProcessUndefined, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CalculateIncrementalMotionProcess(r_model_part, Parameters(R"({"variable_name": "NOTHING"})")),
        "Invalid variable name: NOTHING. Expected DISPLACEMENT or ROTATION.");
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoCalculateIncrementalMotionProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("dummy", 2);
    const CalculateIncrementalMotionProcess process(r_model_part, {R"({"variable_name": "DISPLACEMENT"})"});

    KRATOS_EXPECT_EQ(process.Info(), "CalculateIncrementalMotionProcess");
}

} // namespace Kratos::Testing
