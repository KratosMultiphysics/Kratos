//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "utilities/normal_calculation_utils.h"

namespace Kratos::Testing
{

ModelPart& CreateModelPart(Model& rModel)
{
    PointerVector<Node> node_yz;
    PointerVector<Node> node_xz;
    PointerVector<Node> node_xy;

    GlobalPointersVector<Condition> neighbour_conditions;

    auto& r_model_part = rModel.CreateModelPart("CornerModelPart");
    
    r_model_part.GetProcessInfo()[DOMAIN_SIZE] = 3;

    r_model_part.AddNodalSolutionStepVariable(NORMAL);
    r_model_part.AddNodalSolutionStepVariable(NODAL_PAUX);
    r_model_part.AddNodalSolutionStepVariable(IS_SLIP);

    auto p_properties = r_model_part.CreateNewProperties(0);

    const Node::Pointer p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    const Node::Pointer p_node_2 = r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0);
    const Node::Pointer p_node_3 = r_model_part.CreateNewNode(3, 0.0, 0.0, 1.0);
    const Node::Pointer p_node_4 = r_model_part.CreateNewNode(4, 1.0, 0.0, 0.0);
    const Node::Pointer p_node_5 = r_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
    const Node::Pointer p_node_6 = r_model_part.CreateNewNode(6, 0.0, 1.0, 0.0);

    node_yz.push_back(p_node_1); node_yz.push_back(p_node_2); node_yz.push_back(p_node_3);
    node_xz.push_back(p_node_1); node_xz.push_back(p_node_5); node_xz.push_back(p_node_4);
    node_xy.push_back(p_node_1); node_xy.push_back(p_node_4); node_xy.push_back(p_node_6);

    auto p_condition_yz = r_model_part.CreateNewCondition("SurfaceCondition3D3N", 100, node_yz, p_properties);
    auto p_condition_xz = r_model_part.CreateNewCondition("SurfaceCondition3D3N", 101, node_xz, p_properties);
    auto p_condition_xy = r_model_part.CreateNewCondition("SurfaceCondition3D3N", 102, node_xy, p_properties);

    p_condition_yz->SetValue(PRESSURE, 1.0);
    p_condition_xz->SetValue(PRESSURE, 1.0);
    p_condition_xy->SetValue(PRESSURE, 1.0);

    neighbour_conditions.push_back(p_condition_yz);
    neighbour_conditions.push_back(p_condition_xz);
    neighbour_conditions.push_back(p_condition_xy);

    p_node_1->SetValue(NEIGHBOUR_CONDITIONS, neighbour_conditions);

    return rModel.GetModelPart("CornerModelPart");
}

KRATOS_TEST_CASE_IN_SUITE(NormalCalculationCornerDetection, KratosCoreFastSuite)
{
    Model current_model;
    auto& r_model_part = CreateModelPart(current_model);

    NormalCalculationUtils normal_calculation_utils;
    normal_calculation_utils.CalculateOnSimplex(r_model_part, 3, PRESSURE, 0.0, 10.0);

    const auto& r_corner_node = r_model_part.GetNode(1);
    const array_1d<double, 3> expected_normal(3, 1.0 / 6.0);

    KRATOS_EXPECT_VECTOR_NEAR(r_corner_node.FastGetSolutionStepValue(NORMAL), expected_normal, 1e-10);
    KRATOS_EXPECT_NEAR(r_corner_node.FastGetSolutionStepValue(NODAL_PAUX), 0.5, 1e-10);
    KRATOS_EXPECT_NEAR(r_corner_node.FastGetSolutionStepValue(IS_SLIP), 30.0, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(NormalCalculationCornerDetection_LowMemory, KratosCoreFastSuite)
{
    Model current_model;
    auto& r_model_part = CreateModelPart(current_model);

    NormalCalculationUtils normal_calculation_utils;
    normal_calculation_utils.CalculateOnSimplexLowMemory(r_model_part, 3, PRESSURE, 0.0, 10.0);

    const auto& r_corner_node = r_model_part.GetNode(1);
    const array_1d<double, 3> expected_normal(3, 1.0 / 6.0);
    
    KRATOS_EXPECT_VECTOR_NEAR(r_corner_node.FastGetSolutionStepValue(NORMAL), expected_normal, 1e-10);
    KRATOS_EXPECT_NEAR(r_corner_node.FastGetSolutionStepValue(NODAL_PAUX), 0.5, 1e-10);
    KRATOS_EXPECT_NEAR(r_corner_node.FastGetSolutionStepValue(IS_SLIP), 30.0, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(NormalCalculationCornerDetection_LowMemory_KeepNonHistorical, KratosCoreFastSuite)
{
    Model current_model;
    auto& r_model_part = CreateModelPart(current_model);

    const array_1d<double, 3> test_normal({7.0, -3.5, 1.25});
    double test_nodal_paux = 10.0;

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.SetValue(NORMAL, test_normal);
        r_node.SetValue(NODAL_PAUX, test_nodal_paux);
    }

    NormalCalculationUtils normal_calculation_utils;
    normal_calculation_utils.CalculateOnSimplexLowMemory(r_model_part, 3, PRESSURE, 0.0, 10.0);

    // Test that previous NORMAL and NODAL_PAUX values stored in the non historical database are preserved.
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_VECTOR_NEAR(r_node.GetValue(NORMAL), test_normal, 1e-12);
        KRATOS_EXPECT_NEAR(r_node.GetValue(NODAL_PAUX), test_nodal_paux, 1e-12);
    }
}

} // namespace Kratos::Testing
