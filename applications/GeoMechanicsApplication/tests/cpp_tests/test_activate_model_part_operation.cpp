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
#include "testing/testing.h"

// Application includes
#include "custom_operations/activate_model_part_operation.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(ActivateModelPartOperation, KratosGeoMechanicsFastSuite)
{
    // Create the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");

    // Set up the test model part mesh
    auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({
        "number_of_divisions": 2,
        "element_name": "Element2D3N",
        "condition_name": "LineCondition",
        "create_skin_sub_model_part": true
    })");
    StructuredMeshGeneratorProcess(domain_geometry, r_test_model_part, mesher_parameters).Execute();

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
    Parameters operation_settings(R"({
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
