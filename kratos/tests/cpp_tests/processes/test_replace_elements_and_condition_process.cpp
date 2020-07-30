//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"

#include "elements/distance_calculation_element_simplex.h"
/* Processes */
#include "processes/replace_elements_and_condition_process.h"
#include "utilities/compare_elements_and_conditions_utility.h"

namespace Kratos
{
namespace Testing
{
typedef Node<3> NodeType;

/**
* Checks if the replacement works with triangles (aka 2D geometries)
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess1, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");

    Properties::Pointer p_elem_prop = this_model_part.CreateNewProperties(0);

    // First we create the nodes
    NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    // Now we create the "conditions"
    std::vector<NodeType::Pointer> element_nodes_0 (3);
    element_nodes_0[0] = p_node_1;
    element_nodes_0[1] = p_node_2;
    element_nodes_0[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_0 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_0} );

    std::vector<NodeType::Pointer> element_nodes_1 (3);
    element_nodes_1[0] = p_node_1;
    element_nodes_1[1] = p_node_3;
    element_nodes_1[2] = p_node_4;
    Triangle2D3 <NodeType>::Pointer p_triangle_1 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_1} );

    std::vector<NodeType::Pointer> element_nodes_2 (3);
    element_nodes_2[0] = p_node_2;
    element_nodes_2[1] = p_node_5;
    element_nodes_2[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_2 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_2} );

    std::vector<NodeType::Pointer> element_nodes_3 (3);
    element_nodes_3[0] = p_node_5;
    element_nodes_3[1] = p_node_6;
    element_nodes_3[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_3 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_3} );

    Element::Pointer p_elem_0 = Kratos::make_intrusive<Element>(1, p_triangle_0, p_elem_prop);
    Element::Pointer p_elem_1 = Kratos::make_intrusive<Element>(2, p_triangle_1, p_elem_prop);
    Element::Pointer p_elem_2 = Kratos::make_intrusive<Element>(3, p_triangle_2, p_elem_prop);
    Element::Pointer p_elem_3 = Kratos::make_intrusive<Element>(4, p_triangle_3, p_elem_prop);
    this_model_part.AddElement(p_elem_0);
    this_model_part.AddElement(p_elem_1);
    this_model_part.AddElement(p_elem_2);
    this_model_part.AddElement(p_elem_3);

    // Compute process
    Parameters settings( R"(
    {
        "element_name":"Element2D3N",
        "condition_name": "LineCondition2D2N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(this_model_part, settings);
    process.Execute();

    // Same element type, same geometry type
    std::string component_name;
    for (auto& r_element : this_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_CHECK_EQUAL(component_name, "Element2D3N");
    }
}

/**
* Checks if the replacement works with tetras (aka 3D geometries)
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess2, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");

    Properties::Pointer p_elem_prop = this_model_part.CreateNewProperties(0);

    // First we create the nodes
    NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

    NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_9 = this_model_part.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_10 = this_model_part.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_11 = this_model_part.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_12 = this_model_part.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    // Now we create the "conditions"
    std::vector<NodeType::Pointer> element_nodes_0 (4);
    element_nodes_0[0] = p_node_12;
    element_nodes_0[1] = p_node_10;
    element_nodes_0[2] = p_node_8;
    element_nodes_0[3] = p_node_9;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_0 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_0} );

    std::vector<NodeType::Pointer> element_nodes_1 (4);
    element_nodes_1[0] = p_node_4;
    element_nodes_1[1] = p_node_6;
    element_nodes_1[2] = p_node_9;
    element_nodes_1[3] = p_node_7;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_1 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_1} );

    std::vector<NodeType::Pointer> element_nodes_2 (4);
    element_nodes_2[0] = p_node_11;
    element_nodes_2[1] = p_node_7;
    element_nodes_2[2] = p_node_9;
    element_nodes_2[3] = p_node_8;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_2 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_2} );

    std::vector<NodeType::Pointer> element_nodes_3 (4);
    element_nodes_3[0] = p_node_5;
    element_nodes_3[1] = p_node_3;
    element_nodes_3[2] = p_node_8;
    element_nodes_3[3] = p_node_6;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_3 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_3} );

    std::vector<NodeType::Pointer> element_nodes_4 (4);
    element_nodes_4[0] = p_node_4;
    element_nodes_4[1] = p_node_6;
    element_nodes_4[2] = p_node_7;
    element_nodes_4[3] = p_node_3;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_4 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_4} );

    std::vector<NodeType::Pointer> element_nodes_5 (4);
    element_nodes_5[0] = p_node_2;
    element_nodes_5[1] = p_node_3;
    element_nodes_5[2] = p_node_5;
    element_nodes_5[3] = p_node_6;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_5 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_5} );

    std::vector<NodeType::Pointer> element_nodes_6 (4);
    element_nodes_6[0] = p_node_10;
    element_nodes_6[1] = p_node_9;
    element_nodes_6[2] = p_node_6;
    element_nodes_6[3] = p_node_8;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_6 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_6} );

    std::vector<NodeType::Pointer> element_nodes_7 (4);
    element_nodes_7[0] = p_node_7;
    element_nodes_7[1] = p_node_8;
    element_nodes_7[2] = p_node_3;
    element_nodes_7[3] = p_node_6;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_7 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_7} );

    std::vector<NodeType::Pointer> element_nodes_8 (4);
    element_nodes_8[0] = p_node_7;
    element_nodes_8[1] = p_node_8;
    element_nodes_8[2] = p_node_6;
    element_nodes_8[3] = p_node_9;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_8 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_8} );

    std::vector<NodeType::Pointer> element_nodes_9 (4);
    element_nodes_9[0] = p_node_4;
    element_nodes_9[1] = p_node_1;
    element_nodes_9[2] = p_node_6;
    element_nodes_9[3] = p_node_3;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_9 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_9} );

    std::vector<NodeType::Pointer> element_nodes_10 (4);
    element_nodes_10[0] = p_node_9;
    element_nodes_10[1] = p_node_12;
    element_nodes_10[2] = p_node_11;
    element_nodes_10[3] = p_node_8;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_10 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_10} );

    std::vector<NodeType::Pointer> element_nodes_11 (4);
    element_nodes_11[0] = p_node_3;
    element_nodes_11[1] = p_node_2;
    element_nodes_11[2] = p_node_1;
    element_nodes_11[3] = p_node_6;
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_11 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_11} );

    Element::Pointer p_elem_0 = Kratos::make_intrusive<Element>(1, p_tetrahedra_0, p_elem_prop);
    Element::Pointer p_elem_1 = Kratos::make_intrusive<Element>(2, p_tetrahedra_1, p_elem_prop);
    Element::Pointer p_elem_2 = Kratos::make_intrusive<Element>(3, p_tetrahedra_2, p_elem_prop);
    Element::Pointer p_elem_3 = Kratos::make_intrusive<Element>(4, p_tetrahedra_3, p_elem_prop);
    Element::Pointer p_elem_4 = Kratos::make_intrusive<Element>(5, p_tetrahedra_4, p_elem_prop);
    Element::Pointer p_elem_5 = Kratos::make_intrusive<Element>(6, p_tetrahedra_5, p_elem_prop);
    Element::Pointer p_elem_6 = Kratos::make_intrusive<Element>(7, p_tetrahedra_6, p_elem_prop);
    Element::Pointer p_elem_7 = Kratos::make_intrusive<Element>(8, p_tetrahedra_7, p_elem_prop);
    Element::Pointer p_elem_8 = Kratos::make_intrusive<Element>(9, p_tetrahedra_8, p_elem_prop);
    Element::Pointer p_elem_9 = Kratos::make_intrusive<Element>(10, p_tetrahedra_9, p_elem_prop);
    Element::Pointer p_elem_10 = Kratos::make_intrusive<Element>(11, p_tetrahedra_10, p_elem_prop);
    Element::Pointer p_elem_11 = Kratos::make_intrusive<Element>(12, p_tetrahedra_11, p_elem_prop);
    this_model_part.AddElement(p_elem_0);
    this_model_part.AddElement(p_elem_1);
    this_model_part.AddElement(p_elem_2);
    this_model_part.AddElement(p_elem_3);
    this_model_part.AddElement(p_elem_4);
    this_model_part.AddElement(p_elem_5);
    this_model_part.AddElement(p_elem_6);
    this_model_part.AddElement(p_elem_7);
    this_model_part.AddElement(p_elem_8);
    this_model_part.AddElement(p_elem_9);
    this_model_part.AddElement(p_elem_10);
    this_model_part.AddElement(p_elem_11);

    // Compute process
    Parameters settings( R"(
    {
        "element_name":"Element3D4N",
        "condition_name": "SurfaceCondition3D3N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(this_model_part, settings);
    process.Execute();

    // Same element type
    std::string component_name;
    for (auto& r_element : this_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_CHECK_EQUAL(component_name, "Element3D4N");
    }
}

/**
* Checks if the replacement works when only elements are replaced (not conditions)
*/
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess3, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& this_model_part = current_model.CreateModelPart("Main");

    Properties::Pointer p_elem_prop = this_model_part.CreateNewProperties(0);

    // First we create the nodes
    NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    // Now we create the "conditions"
    std::vector<NodeType::Pointer> element_nodes_0 (3);
    element_nodes_0[0] = p_node_1;
    element_nodes_0[1] = p_node_2;
    element_nodes_0[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_0 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_0} );

    std::vector<NodeType::Pointer> element_nodes_1 (3);
    element_nodes_1[0] = p_node_1;
    element_nodes_1[1] = p_node_3;
    element_nodes_1[2] = p_node_4;
    Triangle2D3 <NodeType>::Pointer p_triangle_1 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_1} );

    std::vector<NodeType::Pointer> element_nodes_2 (3);
    element_nodes_2[0] = p_node_2;
    element_nodes_2[1] = p_node_5;
    element_nodes_2[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_2 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_2} );

    std::vector<NodeType::Pointer> element_nodes_3 (3);
    element_nodes_3[0] = p_node_5;
    element_nodes_3[1] = p_node_6;
    element_nodes_3[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_3 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_3} );

    Element::Pointer p_elem_0 = Kratos::make_intrusive<Element>(1, p_triangle_0, p_elem_prop);
    Element::Pointer p_elem_1 = Kratos::make_intrusive<Element>(2, p_triangle_1, p_elem_prop);
    Element::Pointer p_elem_2 = Kratos::make_intrusive<Element>(3, p_triangle_2, p_elem_prop);
    Element::Pointer p_elem_3 = Kratos::make_intrusive<Element>(4, p_triangle_3, p_elem_prop);
    this_model_part.AddElement(p_elem_0);
    this_model_part.AddElement(p_elem_1);
    this_model_part.AddElement(p_elem_2);
    this_model_part.AddElement(p_elem_3);

    Condition::Pointer p_cond_0 = Kratos::make_intrusive<Condition>(1, p_triangle_0, p_elem_prop);
    Condition::Pointer p_cond_1 = Kratos::make_intrusive<Condition>(2, p_triangle_1, p_elem_prop);
    Condition::Pointer p_cond_2 = Kratos::make_intrusive<Condition>(3, p_triangle_2, p_elem_prop);
    Condition::Pointer p_cond_3 = Kratos::make_intrusive<Condition>(4, p_triangle_3, p_elem_prop);
    this_model_part.AddCondition(p_cond_0);
    this_model_part.AddCondition(p_cond_1);
    this_model_part.AddCondition(p_cond_2);
    this_model_part.AddCondition(p_cond_3);

    // Compute process
    Parameters settings( R"(
    {
        "element_name":"Element2D3N",
        "condition_name": ""
    }  )" );

    ReplaceElementsAndConditionsProcess process(this_model_part, settings);
    process.Execute();

    // Same element type, same geometry type
    std::string component_name;
    for (auto& r_element : this_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_CHECK_EQUAL(component_name, "Element2D3N");
    }

    // check that conditions were NOT replaced
    KRATOS_CHECK_EQUAL(this_model_part.pGetCondition(1), p_cond_0);
    KRATOS_CHECK_EQUAL(this_model_part.pGetCondition(2), p_cond_1);
    KRATOS_CHECK_EQUAL(this_model_part.pGetCondition(3), p_cond_2);
    KRATOS_CHECK_EQUAL(this_model_part.pGetCondition(4), p_cond_3);
}

// We check that if a submodelpart is given only specified elements are replaced in the rest
// of model parts
KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess4, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& parent_model_part = current_model.CreateModelPart("Parent");
    ModelPart& this_model_part = parent_model_part.CreateSubModelPart("Main");
    ModelPart& sister_model_part = parent_model_part.CreateSubModelPart("Sister");

    Properties::Pointer p_elem_prop = this_model_part.CreateNewProperties(0);

    // First we create the nodes
    NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    // Now we create the "geometries"
    std::vector<NodeType::Pointer> element_nodes_0 (3);
    element_nodes_0[0] = p_node_1;
    element_nodes_0[1] = p_node_2;
    element_nodes_0[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_0 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_0} );

    std::vector<NodeType::Pointer> element_nodes_1 (3);
    element_nodes_1[0] = p_node_1;
    element_nodes_1[1] = p_node_3;
    element_nodes_1[2] = p_node_4;
    Triangle2D3 <NodeType>::Pointer p_triangle_1 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_1} );

    std::vector<NodeType::Pointer> element_nodes_2 (3);
    element_nodes_2[0] = p_node_2;
    element_nodes_2[1] = p_node_5;
    element_nodes_2[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_2 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_2} );

    std::vector<NodeType::Pointer> element_nodes_3 (3);
    element_nodes_3[0] = p_node_5;
    element_nodes_3[1] = p_node_6;
    element_nodes_3[2] = p_node_3;
    Triangle2D3 <NodeType>::Pointer p_triangle_3 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_3} );

    // Elements are created
    Element::Pointer p_elem_0 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(1, p_triangle_0, p_elem_prop);
    Element::Pointer p_elem_1 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(2, p_triangle_1, p_elem_prop);
    Element::Pointer p_elem_2 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(3, p_triangle_2, p_elem_prop);
    Element::Pointer p_elem_3 = Kratos::make_intrusive<DistanceCalculationElementSimplex<2>>(4, p_triangle_3, p_elem_prop);

    // Element with id 4 is present in this model part and should therefore not be replaced
    parent_model_part.AddElement(p_elem_0);
    parent_model_part.AddElement(p_elem_1);
    parent_model_part.AddElement(p_elem_2);
    parent_model_part.AddElement(p_elem_3);

    this_model_part.AddElement(p_elem_0);
    this_model_part.AddElement(p_elem_1);
    this_model_part.AddElement(p_elem_2);

    sister_model_part.AddElement(p_elem_0);
    sister_model_part.AddElement(p_elem_1);
    sister_model_part.AddElement(p_elem_3);

    // Compute process
    Parameters settings( R"(
    {
        "element_name":"Element2D3N",
        "condition_name": "LineCondition2D2N"
    }  )" );

    ReplaceElementsAndConditionsProcess process(this_model_part, settings);
    process.Execute();

    std::vector<int> modified_elems_ids = {1 ,2 ,3};
    std::string component_name;
    for (auto& r_element : this_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        KRATOS_CHECK_EQUAL(component_name, "Element2D3N");
    }

    for (auto& r_element : sister_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        if (std::find(modified_elems_ids.begin(), modified_elems_ids.end(), r_element.Id()) != modified_elems_ids.end()) {
            KRATOS_CHECK_EQUAL(component_name, "Element2D3N");
        }
        else {
            KRATOS_CHECK_EQUAL(component_name, "DistanceCalculationElementSimplex2D3N");
        }
        
    }

    for (auto& r_element : parent_model_part.Elements()) {
        CompareElementsAndConditionsUtility::GetRegisteredName(r_element, component_name);
        if (std::find(modified_elems_ids.begin(), modified_elems_ids.end(), r_element.Id()) != modified_elems_ids.end()) {
            KRATOS_CHECK_EQUAL(component_name, "Element2D3N");
        }
        else {
            KRATOS_CHECK_EQUAL(component_name, "DistanceCalculationElementSimplex2D3N");
        }
        
    }
}

} // namespace Testing
}  // namespace Kratos.
