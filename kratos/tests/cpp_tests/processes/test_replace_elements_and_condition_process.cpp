//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "elements/distance_calculation_element_simplex.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/replace_elements_and_condition_process.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos::Testing
{

void GenerateMixedModelPart(ModelPart& rModelPart)
{
    Properties::Pointer p_prop = rModelPart.CreateNewProperties(0);

    // First we create the nodes
    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    auto p_node_4 = rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    auto p_node_5 = rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    auto p_node_6 = rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    // Now we create the "elements"
    {
        std::vector<Node::Pointer> element_nodes_0 ({p_node_1, p_node_2, p_node_3});
        Triangle2D3 <Node>::Pointer p_geom_0 = Kratos::make_shared<Triangle2D3 <Node>>( PointerVector<Node>{element_nodes_0} );

        std::vector<Node::Pointer> element_nodes_1 ({p_node_1, p_node_3, p_node_4});
        Triangle2D3 <Node>::Pointer p_geom_1 = Kratos::make_shared<Triangle2D3 <Node>>( PointerVector<Node>{element_nodes_1} );

        std::vector<Node::Pointer> element_nodes_2 ({p_node_2, p_node_5, p_node_6, p_node_3});
        Quadrilateral2D4 <Node>::Pointer p_geom_2 = Kratos::make_shared<Quadrilateral2D4 <Node>>( PointerVector<Node>{element_nodes_2} );

        auto p_elem_0 = Kratos::make_intrusive<Element>(1, p_geom_0, p_prop);
        auto p_elem_1 = Kratos::make_intrusive<Element>(2, p_geom_1, p_prop);
        auto p_elem_2 = Kratos::make_intrusive<Element>(3, p_geom_2, p_prop);
        rModelPart.AddElement(p_elem_0);
        rModelPart.AddElement(p_elem_1);
        rModelPart.AddElement(p_elem_2);
    }

    // Now we create the "conditions"
    {
        std::vector<Node::Pointer> conditions_nodes_0 ({p_node_1, p_node_2});
        Line2D2 <Node>::Pointer p_geom_0 = Kratos::make_shared<Line2D2 <Node>>( PointerVector<Node>{conditions_nodes_0} );

        std::vector<Node::Pointer> conditions_nodes_1 ({p_node_1, p_node_5});
        Line2D2 <Node>::Pointer p_geom_1 = Kratos::make_shared<Line2D2 <Node>>( PointerVector<Node>{conditions_nodes_1} );

        std::vector<Node::Pointer> conditions_nodes_2 ({p_node_5, p_node_6});
        Line2D2 <Node>::Pointer p_geom_2 = Kratos::make_shared<Line2D2 <Node>>( PointerVector<Node>{conditions_nodes_2} );

        std::vector<Node::Pointer> conditions_nodes_3 ({p_node_6, p_node_3});
        Line2D2 <Node>::Pointer p_geom_3 = Kratos::make_shared<Line2D2 <Node>>( PointerVector<Node>{conditions_nodes_3} );

        std::vector<Node::Pointer> conditions_nodes_4 ({p_node_3, p_node_4});
        Line2D2 <Node>::Pointer p_geom_4 = Kratos::make_shared<Line2D2 <Node>>( PointerVector<Node>{conditions_nodes_4} );

        std::vector<Node::Pointer> conditions_nodes_5 ({p_node_4, p_node_1});
        Line2D2 <Node>::Pointer p_geom_5 = Kratos::make_shared<Line2D2 <Node>>( PointerVector<Node>{conditions_nodes_5} );

        auto p_cond_0 = Kratos::make_intrusive<Condition>(1, p_geom_0, p_prop);
        auto p_cond_1 = Kratos::make_intrusive<Condition>(2, p_geom_1, p_prop);
        auto p_cond_2 = Kratos::make_intrusive<Condition>(3, p_geom_2, p_prop);
        auto p_cond_3 = Kratos::make_intrusive<Condition>(4, p_geom_3, p_prop);
        auto p_cond_4 = Kratos::make_intrusive<Condition>(5, p_geom_4, p_prop);
        auto p_cond_5 = Kratos::make_intrusive<Condition>(6, p_geom_5, p_prop);
        rModelPart.AddCondition(p_cond_0);
        rModelPart.AddCondition(p_cond_1);
        rModelPart.AddCondition(p_cond_2);
        rModelPart.AddCondition(p_cond_3);
        rModelPart.AddCondition(p_cond_4);
        rModelPart.AddCondition(p_cond_5);
    }
}

/**
* Checks if the replacement works with triangles (aka 2D geometries)
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess1, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // Generate a simple mesh
    CppTestsUtilities::CreateTestModelPartTriangle2D3N(r_model_part);

    const double pressure = 3.0;
    for (auto& r_element : r_model_part.Elements()) {
        r_element.SetValue(PRESSURE, pressure);
        r_element.Set(VISITED);
    }

    // Compute process
    Parameters settings( R"(
    {
        "element_name"   : "Element2D3N",
        "condition_name" : "LineCondition2D2N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(r_model_part, settings);
    process.Execute();

    // Same element type, same geometry type
    std::string component_name;
    for (auto& r_element : r_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_EXPECT_EQ(component_name, "Element2D3N");
        KRATOS_EXPECT_NEAR(r_element.GetValue(PRESSURE), pressure, 1e-6);
        KRATOS_EXPECT_TRUE(r_element.Is(VISITED));
    }
}

/**
* Checks if the replacement works with tetras (aka 3D geometries)
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess2, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // Generate a simple mesh
    CppTestsUtilities::CreateTestModelPartTetrahedra3D4N(r_model_part);

    // Compute process
    Parameters settings( R"(
    {
        "element_name"   : "Element3D4N",
        "condition_name" : "SurfaceCondition3D3N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(r_model_part, settings);
    process.Execute();

    // Same element type
    std::string component_name;
    for (auto& r_element : r_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_EXPECT_EQ(component_name, "Element3D4N");
    }
}

/**
* Checks if the replacement works when only elements are replaced (not conditions)
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess3, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // Generate a simple mesh
    CppTestsUtilities::CreateTestModelPartTriangle2D3N(r_model_part);

    auto p_prop = r_model_part.pGetProperties(0);
    auto p_cond_0 = Kratos::make_intrusive<Condition>(1, r_model_part.pGetElement(1)->pGetGeometry(), p_prop);
    auto p_cond_1 = Kratos::make_intrusive<Condition>(2, r_model_part.pGetElement(2)->pGetGeometry(), p_prop);
    auto p_cond_2 = Kratos::make_intrusive<Condition>(3, r_model_part.pGetElement(3)->pGetGeometry(), p_prop);
    auto p_cond_3 = Kratos::make_intrusive<Condition>(4, r_model_part.pGetElement(4)->pGetGeometry(), p_prop);
    r_model_part.AddCondition(p_cond_0);
    r_model_part.AddCondition(p_cond_1);
    r_model_part.AddCondition(p_cond_2);
    r_model_part.AddCondition(p_cond_3);

    // Compute process
    Parameters settings( R"(
    {
        "element_name"   : "Element2D3N",
        "condition_name" : ""
    }  )" );

    ReplaceElementsAndConditionsProcess process(r_model_part, settings);
    process.Execute();

    // Same element type, same geometry type
    std::string component_name;
    for (auto& r_element : r_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_EXPECT_EQ(component_name, "Element2D3N");
    }

    // Check that conditions were NOT replaced
    KRATOS_EXPECT_EQ(r_model_part.pGetCondition(1), p_cond_0);
    KRATOS_EXPECT_EQ(r_model_part.pGetCondition(2), p_cond_1);
    KRATOS_EXPECT_EQ(r_model_part.pGetCondition(3), p_cond_2);
    KRATOS_EXPECT_EQ(r_model_part.pGetCondition(4), p_cond_3);
}

// We check that if a submodelpart is given only specified elements are replaced in the rest
// of model parts
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess4, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_parent_model_part = current_model.CreateModelPart("Parent");
    ModelPart& r_model_part = r_parent_model_part.CreateSubModelPart("Main");
    ModelPart& r_sister_model_part = r_parent_model_part.CreateSubModelPart("Sister");

    using Node = Node;

    // Generate a simple mesh
    Properties::Pointer p_prop = r_model_part.CreateNewProperties(0);

    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    auto p_node_5 = r_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    auto p_node_6 = r_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    // Now we create the "conditions"
    std::vector<Node::Pointer> element_nodes_0 ({p_node_1, p_node_2, p_node_3});
    Triangle2D3 <Node>::Pointer p_geom_1 = Kratos::make_shared<Triangle2D3 <Node>>( PointerVector<Node>{element_nodes_0} );

    std::vector<Node::Pointer> element_nodes_1 ({p_node_1, p_node_3, p_node_4});
    Triangle2D3 <Node>::Pointer p_geom_2 = Kratos::make_shared<Triangle2D3 <Node>>( PointerVector<Node>{element_nodes_1} );

    std::vector<Node::Pointer> element_nodes_2 ({p_node_2, p_node_5, p_node_3});
    Triangle2D3 <Node>::Pointer p_geom_3 = Kratos::make_shared<Triangle2D3 <Node>>( PointerVector<Node>{element_nodes_2} );

    std::vector<Node::Pointer> element_nodes_3 ({p_node_5, p_node_6, p_node_3});
    Triangle2D3 <Node>::Pointer p_geom_4 = Kratos::make_shared<Triangle2D3 <Node>>( PointerVector<Node>{element_nodes_3} );

    // Elements are created
    auto p_elem_0 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(1, p_geom_1, p_prop);
    auto p_elem_1 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(2, p_geom_1, p_prop);
    auto p_elem_2 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(3, p_geom_1, p_prop);
    auto p_elem_3 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(4, p_geom_1, p_prop);

    // Element with id 4 is present in this model part and should therefore not be replaced
    r_parent_model_part.AddElement(p_elem_0);
    r_parent_model_part.AddElement(p_elem_1);
    r_parent_model_part.AddElement(p_elem_2);
    r_parent_model_part.AddElement(p_elem_3);

    r_model_part.AddElement(p_elem_0);
    r_model_part.AddElement(p_elem_1);
    r_model_part.AddElement(p_elem_2);

    r_sister_model_part.AddElement(p_elem_0);
    r_sister_model_part.AddElement(p_elem_1);
    r_sister_model_part.AddElement(p_elem_3);

    // Compute process
    Parameters settings( R"(
    {
        "element_name"  : "Element2D3N",
        "condition_name": "LineCondition2D2N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(r_model_part, settings);
    process.Execute();

    std::vector<int> modified_elems_ids = {1 ,2 ,3};
    std::string component_name;
    for (auto& r_element : r_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_EXPECT_EQ(component_name, "Element2D3N");
    }

    for (auto& r_element : r_sister_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        if (std::find(modified_elems_ids.begin(), modified_elems_ids.end(), r_element.Id()) != modified_elems_ids.end()) {
            KRATOS_EXPECT_EQ(component_name, "Element2D3N");
        } else {
            KRATOS_EXPECT_EQ(component_name, "DistanceCalculationElementSimplex2D3N");
        }
    }

    for (auto& r_element : r_parent_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        if (std::find(modified_elems_ids.begin(), modified_elems_ids.end(), r_element.Id()) != modified_elems_ids.end()) {
            KRATOS_EXPECT_EQ(component_name, "Element2D3N");
        }else {
            KRATOS_EXPECT_EQ(component_name, "DistanceCalculationElementSimplex2D3N");
        }
    }
}

/**
* Checks if the replacement works with triangles combined geometries
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess5, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // Generated mixed model part
    GenerateMixedModelPart(r_model_part);

    // Compute process
    Parameters settings( R"(
    {
        "element_name"   : "Element2D3N;Element2D4N",
        "condition_name" : "LineCondition2D2N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(r_model_part, settings);
    process.Execute();

    // Same element type, same geometry type
    std::string component_name;
    for (auto& r_element : r_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        if (r_element.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
            KRATOS_EXPECT_EQ(component_name, "Element2D3N");
        } else if (r_element.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4) {
            KRATOS_EXPECT_EQ(component_name, "Element2D4N");
        }
    }
    for (auto& r_condition : r_model_part.Conditions()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_condition, component_name);
        KRATOS_EXPECT_EQ(component_name, "LineCondition2D2N");
    }
}

/**
* Checks if the replacement works with triangles combined geometries
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess6, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    // Generated mixed model part
    GenerateMixedModelPart(r_model_part);

    // Compute process
    Parameters settings( R"(
    {
        "element_name"   : "Element#D#N",
        "condition_name" : "LineCondition2D2N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(r_model_part, settings);
    process.Execute();

    // Same element type, same geometry type
    std::string component_name;
    for (auto& r_element : r_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        if (r_element.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
            KRATOS_EXPECT_EQ(component_name, "Element2D3N");
        } else if (r_element.GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4) {
            KRATOS_EXPECT_EQ(component_name, "Element2D4N");
        }
    }
    for (auto& r_condition : r_model_part.Conditions()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_condition, component_name);
        KRATOS_EXPECT_EQ(component_name, "LineCondition2D2N");
    }
}

} // namespace Kratos::Testing
