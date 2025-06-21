// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// Project includes
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/expect.h"
#include "processes/structured_mesh_generator_process.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

// Application includes
#include "custom_operations/activate_model_part_operation.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ActivateModelPartOperation, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Create the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    auto  p_properties      = r_test_model_part.CreateNewProperties(0);

    const auto p_node1 = r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    const auto p_node2 = r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    const auto p_node3 = r_test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    const auto p_node4 = r_test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);

    r_test_model_part.AddElement(make_intrusive<Element>(
        1, std::make_shared<Triangle2D3<Node>>(p_node1, p_node2, p_node3), p_properties));
    r_test_model_part.AddElement(make_intrusive<Element>(
        3, std::make_shared<Triangle2D3<Node>>(p_node1, p_node3, p_node4), p_properties));

    // Deactivate all the model part entities
    for (auto& r_node : r_test_model_part.Nodes()) {
        r_node.Set(ACTIVE, false);
    }
    for (auto& r_element : r_test_model_part.Elements()) {
        r_element.Set(ACTIVE, false);
    }
    for (auto& r_condition : r_test_model_part.Conditions()) {
        r_condition.Set(ACTIVE, false);
    }

    // Create and execute the tested operation
    Parameters                 operation_settings(R"({
        "model_part_name" : "TestModelPart"
    })");
    ActivateModelPartOperation test_operation(test_model, operation_settings);
    test_operation.Execute();

    // Check that all model part entities are now active
    for (const auto& r_node : r_test_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.Is(ACTIVE))
    }
    for (const auto& r_element : r_test_model_part.Elements()) {
        KRATOS_EXPECT_TRUE(r_element.Is(ACTIVE))
    }
    for (const auto& r_condition : r_test_model_part.Conditions()) {
        KRATOS_EXPECT_TRUE(r_condition.Is(ACTIVE))
    }
}

} // namespace Kratos::Testing
