//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos
{
namespace CppTestsUtilities
{
void Create2DGeometry(
    ModelPart& rModelPart,
    const std::string& rEntityName,
    const bool Initialize,
    const bool Elements
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    if (Elements) {
        rModelPart.CreateNewElement(rEntityName, 1, {{1,2,3}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 2, {{1,3,4}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 3, {{2,5,3}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 4, {{5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_elem : rModelPart.Elements())
                r_elem.Initialize(r_process_info);
        }
    } else {
        rModelPart.CreateNewCondition(rEntityName, 1, {{1,2,3}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 2, {{1,3,4}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 3, {{2,5,3}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 4, {{5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_cond : rModelPart.Conditions())
                r_cond.Initialize(r_process_info);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CreateTestModelPartTriangle2D3N(ModelPart& rModelPart)
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    using NodeType = Node;

    // Clear the model part
    rModelPart.Clear();

    // First we create the nodes
    auto p_node_1 = rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_2 = rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    auto p_node_4 = rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    auto p_node_5 = rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    auto p_node_6 = rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    // Now we create the "elements"
    std::vector<NodeType::Pointer> element_nodes_0 ({p_node_1, p_node_2, p_node_3});
    Triangle2D3 <NodeType>::Pointer p_triangle_0 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_0} );

    std::vector<NodeType::Pointer> element_nodes_1 ({p_node_1, p_node_3, p_node_4});
    Triangle2D3 <NodeType>::Pointer p_triangle_1 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_1} );

    std::vector<NodeType::Pointer> element_nodes_2 ({p_node_2, p_node_5, p_node_3});
    Triangle2D3 <NodeType>::Pointer p_triangle_2 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_2} );

    std::vector<NodeType::Pointer> element_nodes_3 ({p_node_5, p_node_6, p_node_3});
    Triangle2D3 <NodeType>::Pointer p_triangle_3 = Kratos::make_shared<Triangle2D3 <NodeType>>( PointerVector<NodeType>{element_nodes_3} );

    auto p_elem_0 = Kratos::make_intrusive<Element>(1, p_triangle_0, p_prop);
    auto p_elem_1 = Kratos::make_intrusive<Element>(2, p_triangle_1, p_prop);
    auto p_elem_2 = Kratos::make_intrusive<Element>(3, p_triangle_2, p_prop);
    auto p_elem_3 = Kratos::make_intrusive<Element>(4, p_triangle_3, p_prop);
    rModelPart.AddElement(p_elem_0);
    rModelPart.AddElement(p_elem_1);
    rModelPart.AddElement(p_elem_2);
    rModelPart.AddElement(p_elem_3);
}

/***********************************************************************************/
/***********************************************************************************/

void Create2DQuadrilateralsGeometry(
    ModelPart& rModelPart,
    const std::string& rEntityName,
    const bool Initialize,
    const bool Elements
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

    if (Elements) {
        rModelPart.CreateNewElement(rEntityName, 1, {{1,2,3,4}}, p_prop);
        rModelPart.CreateNewElement(rEntityName, 2, {{2,5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_elem : rModelPart.Elements())
                r_elem.Initialize(r_process_info);
        }
    } else {
        rModelPart.CreateNewCondition(rEntityName, 1, {{1,2,3,4}}, p_prop);
        rModelPart.CreateNewCondition(rEntityName, 2, {{2,5,6,3}}, p_prop);

        // Initialize Elements
        if (Initialize) {
            const auto& r_process_info = rModelPart.GetProcessInfo();
            for (auto& r_cond : rModelPart.Conditions())
                r_cond.Initialize(r_process_info);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Create3DGeometry(
    ModelPart& rModelPart,
    const std::string& rElementName,
    const bool Initialize
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

    rModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    rModelPart.CreateNewElement(rElementName, 1, {{12,10,8,9}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 2, {{4,6,9,7}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 3, {{11,7,9,8}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 4, {{5,3,8,6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 5, {{4,6,7,3}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 6, {{2,3,5,6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 7, {{10,9,6,8}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 8, {{7,8,3,6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 9, {{7,8,6,9}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 10, {{4,1,6,3}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 11, {{9,12,11,8}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 12, {{3,2,1,6}}, p_prop);

    // Initialize Elements
    if (Initialize) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& r_elem : rModelPart.Elements())
            r_elem.Initialize(r_process_info);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CreateTestModelPartTetrahedra3D4N(ModelPart& rModelPart)
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    using NodeType = Node;

    // Clear the model part
    rModelPart.Clear();

    // First we create the nodes
    auto p_node_1 = rModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    auto p_node_2 = rModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    auto p_node_3 = rModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    auto p_node_4 = rModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    auto p_node_5 = rModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    auto p_node_6 = rModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

    auto p_node_7 = rModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    auto p_node_8 = rModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    auto p_node_9 = rModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    auto p_node_10 = rModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    auto p_node_11 = rModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    auto p_node_12 = rModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    // Now we create the "elements"
    std::vector<NodeType::Pointer> element_nodes_0 ({p_node_12, p_node_10, p_node_8, p_node_9});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_0 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_0} );

    std::vector<NodeType::Pointer> element_nodes_1 ({p_node_4, p_node_6, p_node_9, p_node_7});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_1 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_1} );

    std::vector<NodeType::Pointer> element_nodes_2 ({p_node_11, p_node_7, p_node_9, p_node_8});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_2 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_2} );

    std::vector<NodeType::Pointer> element_nodes_3 ({p_node_5, p_node_3, p_node_8, p_node_6});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_3 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_3} );

    std::vector<NodeType::Pointer> element_nodes_4 ({p_node_4, p_node_6, p_node_7, p_node_3});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_4 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_4} );

    std::vector<NodeType::Pointer> element_nodes_5 ({p_node_2, p_node_3, p_node_5, p_node_6});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_5 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_5} );

    std::vector<NodeType::Pointer> element_nodes_6 ({p_node_10, p_node_9, p_node_6, p_node_8});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_6 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_6} );

    std::vector<NodeType::Pointer> element_nodes_7 ({p_node_7, p_node_8, p_node_3, p_node_6});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_7 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_7} );

    std::vector<NodeType::Pointer> element_nodes_8 ({p_node_7, p_node_8, p_node_6, p_node_9});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_8 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_8} );

    std::vector<NodeType::Pointer> element_nodes_9 ({p_node_4, p_node_1, p_node_6, p_node_3});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_9 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_9} );

    std::vector<NodeType::Pointer> element_nodes_10 ({p_node_9, p_node_12, p_node_11, p_node_8});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_10 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_10} );

    std::vector<NodeType::Pointer> element_nodes_11 ({p_node_3, p_node_2, p_node_1, p_node_6});
    Tetrahedra3D4 <NodeType>::Pointer p_tetrahedra_11 = Kratos::make_shared<Tetrahedra3D4 <NodeType>>( PointerVector<NodeType>{element_nodes_11} );

    auto p_elem_0 = Kratos::make_intrusive<Element>(1, p_tetrahedra_0, p_prop);
    auto p_elem_1 = Kratos::make_intrusive<Element>(2, p_tetrahedra_1, p_prop);
    auto p_elem_2 = Kratos::make_intrusive<Element>(3, p_tetrahedra_2, p_prop);
    auto p_elem_3 = Kratos::make_intrusive<Element>(4, p_tetrahedra_3, p_prop);
    auto p_elem_4 = Kratos::make_intrusive<Element>(5, p_tetrahedra_4, p_prop);
    auto p_elem_5 = Kratos::make_intrusive<Element>(6, p_tetrahedra_5, p_prop);
    auto p_elem_6 = Kratos::make_intrusive<Element>(7, p_tetrahedra_6, p_prop);
    auto p_elem_7 = Kratos::make_intrusive<Element>(8, p_tetrahedra_7, p_prop);
    auto p_elem_8 = Kratos::make_intrusive<Element>(9, p_tetrahedra_8, p_prop);
    auto p_elem_9 = Kratos::make_intrusive<Element>(10, p_tetrahedra_9, p_prop);
    auto p_elem_10 = Kratos::make_intrusive<Element>(11, p_tetrahedra_10, p_prop);
    auto p_elem_11 = Kratos::make_intrusive<Element>(12, p_tetrahedra_11, p_prop);
    rModelPart.AddElement(p_elem_0);
    rModelPart.AddElement(p_elem_1);
    rModelPart.AddElement(p_elem_2);
    rModelPart.AddElement(p_elem_3);
    rModelPart.AddElement(p_elem_4);
    rModelPart.AddElement(p_elem_5);
    rModelPart.AddElement(p_elem_6);
    rModelPart.AddElement(p_elem_7);
    rModelPart.AddElement(p_elem_8);
    rModelPart.AddElement(p_elem_9);
    rModelPart.AddElement(p_elem_10);
    rModelPart.AddElement(p_elem_11);
}

/***********************************************************************************/
/***********************************************************************************/

void Create3DHexahedraGeometry(
    ModelPart& rModelPart,
    const std::string& rElementName,
    const bool Initialize
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

    rModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    rModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    rModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    rModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    rModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    rModelPart.CreateNewElement(rElementName, 1, {{5,8,6,2,3,7,4,1}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 2, {{8,12,10,6,7,11,9,4}}, p_prop);

    // Initialize Elements
    if (Initialize) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& r_elem : rModelPart.Elements())
            r_elem.Initialize(r_process_info);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void Create3DQuadraticGeometry(
    ModelPart& rModelPart, 
    const std::string& rElementName,
    const bool Initialize
    )
{
    Properties::Pointer p_prop = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);

    // First we create the nodes
    rModelPart.CreateNewNode(1, 1.0, 0.0, 1.0);
    rModelPart.CreateNewNode(2, 0.5, 0.0, 1.0);
    rModelPart.CreateNewNode(3, 1.0, 0.5, 1.0);
    rModelPart.CreateNewNode(4, 1.0, 0.0, 0.5);
    rModelPart.CreateNewNode(5, 0.5, 0.0, 0.5);
    rModelPart.CreateNewNode(6, 1.0, 0.5, 0.5);
    rModelPart.CreateNewNode(7, 0.5, 0.5, 1.0);
    rModelPart.CreateNewNode(8, 0.5, 0.5, 0.5);
    rModelPart.CreateNewNode(9, 0.2019246055, 0.3959160307, 0.6930668948);
    rModelPart.CreateNewNode(10, 0.7019246055, 0.8959160307, 0.6930668948);
    rModelPart.CreateNewNode(11, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(12, 0.0, 0.0, 1.0);
    rModelPart.CreateNewNode(13, 1.0, 1.0, 1.0);
    rModelPart.CreateNewNode(14, 0.5, 0.0, 0.0);
    rModelPart.CreateNewNode(15, 1.0, 0.5, 0.0);
    rModelPart.CreateNewNode(16, 0.5, 1.0, 1.0);
    rModelPart.CreateNewNode(17, 0.0, 0.5, 1.0);
    rModelPart.CreateNewNode(18, 0.0, 0.0, 0.5);
    rModelPart.CreateNewNode(19, 1.0, 1.0, 0.5);
    rModelPart.CreateNewNode(20, 0.4038492111, 0.7918320615, 0.3861337896);
    rModelPart.CreateNewNode(21, 0.2019246055, 0.3959160307, 0.1930668948);
    rModelPart.CreateNewNode(22, 0.5, 0.5, 0.0);
    rModelPart.CreateNewNode(23, 0.5, 1.0, 0.5);
    rModelPart.CreateNewNode(24, 0.0, 0.5, 0.5);
    rModelPart.CreateNewNode(25, 0.2019246055, 0.8959160307, 0.6930668948);
    rModelPart.CreateNewNode(26, 0.7019246055, 0.8959160307, 0.1930668948);
    rModelPart.CreateNewNode(27, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(28, 1.0, 1.0, 0.0);
    rModelPart.CreateNewNode(29, 0.0, 1.0, 1.0);
    rModelPart.CreateNewNode(30, 0.2019246055, 0.8959160307, 0.1930668948);
    rModelPart.CreateNewNode(31, 0.5, 1.0, 0.0);
    rModelPart.CreateNewNode(32, 0.0, 0.5, 0.0);
    rModelPart.CreateNewNode(33, 0.0, 1.0, 0.5);
    rModelPart.CreateNewNode(34, 0.0, 1.0, 0.0);

    // Now we create the elements
    rModelPart.CreateNewElement(rElementName, 1,  {{27, 11, 28, 13, 14, 15, 22,  8,  6, 19}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 2,  {{34, 12, 27, 20, 24, 18, 32, 30,  9, 21}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 3,  {{27, 34, 20, 28, 32, 30, 21, 22, 31, 26}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 4,  {{34, 20, 28, 29, 30, 26, 31, 33, 25, 23}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 5,  {{34, 20, 29, 12, 30, 25, 33, 24,  9, 17}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 6,  {{20, 27, 28, 13, 21, 22, 26, 10,  8, 19}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 7,  {{28, 20, 13, 29, 26, 10, 19, 23, 25, 16}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 8,  {{20, 13, 29, 12, 10, 16, 25,  9,  7, 17}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 9,  {{27,  1, 11, 13,  5,  4, 14,  8,  3,  6}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 10, {{12, 13,  1, 27,  7,  3,  2, 18,  8,  5}}, p_prop);
    rModelPart.CreateNewElement(rElementName, 11, {{12, 27, 20, 13, 18, 21,  9,  7,  8, 10}}, p_prop);

    // Initialize Elements
    if (Initialize) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        for (auto& r_elem : rModelPart.Elements())
            r_elem.Initialize(r_process_info);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CreateSphereTriangularMesh(
    ModelPart& rModelPart,
    const std::string& rConditionName,
    const double Radius,
    const std::array<double, 3>& rCenter
    )
{
    const double scale = Radius/0.5;
    Properties::Pointer p_prop = nullptr;

    // Get the Initial Id
    const std::size_t initial_node_id = block_for_each<MaxReduction<std::size_t>>(rModelPart.GetRootModelPart().Nodes(), [](Node& rNode) {
        return rNode.Id();
    });

    // First we create the nodes
    rModelPart.CreateNewNode(initial_node_id + 1, scale * 0.16372, scale * 0.18066, scale * -0.43653);
    rModelPart.CreateNewNode(initial_node_id + 2, 0.0, scale * 0.15451, scale * -0.47553);
    rModelPart.CreateNewNode(initial_node_id + 3, 0.0, scale * 0.29389, scale * -0.40451);
    rModelPart.CreateNewNode(initial_node_id + 4, scale * -0, scale * -0, scale * -0.5);
    rModelPart.CreateNewNode(initial_node_id + 5, scale * 0.33997, scale * -0, scale * -0.36663);
    rModelPart.CreateNewNode(initial_node_id + 6, scale * 0.23141, scale * 0.37802, scale * -0.23141);
    rModelPart.CreateNewNode(initial_node_id + 7, 0.0, scale * 0.40451, scale * -0.29389);
    rModelPart.CreateNewNode(initial_node_id + 8, scale * -0.16372, scale * 0.18066, scale * -0.43653);
    rModelPart.CreateNewNode(initial_node_id + 9, scale * 0.16372, scale * -0.18066, scale * -0.43653);
    rModelPart.CreateNewNode(initial_node_id + 10, 0.0, scale * -0.15451, scale * -0.47553);
    rModelPart.CreateNewNode(initial_node_id + 11, scale * 0.43653, scale * 0.18066, scale * -0.16372);
    rModelPart.CreateNewNode(initial_node_id + 12, 0.0, scale * 0.47553, scale * -0.15451);
    rModelPart.CreateNewNode(initial_node_id + 13, scale * -0.23141, scale * 0.37802, scale * -0.23141);
    rModelPart.CreateNewNode(initial_node_id + 14, scale * -0.16372, scale * -0.18066, scale * -0.43653);
    rModelPart.CreateNewNode(initial_node_id + 15, 0.0, scale * -0.29389, scale * -0.40451);
    rModelPart.CreateNewNode(initial_node_id + 16, scale * 0.29389, scale * 0.40451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 17, scale * 0.40451, scale * 0.29389, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 18, scale * 0.15451, scale * 0.47553, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 19, scale * 0.43653, scale * -0.18066, scale * -0.16372);
    rModelPart.CreateNewNode(initial_node_id + 20, scale * 0.47553, scale * 0.15451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 21, 0.0, scale * 0.5, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 22, scale * -0.36663, 0.0, scale * -0.33997);
    rModelPart.CreateNewNode(initial_node_id + 23, scale * 0.5, scale * -0, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 24, scale * 0.23141, scale * -0.37802, scale * -0.23141);
    rModelPart.CreateNewNode(initial_node_id + 25, scale * -0.15451, scale * 0.47553, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 26, 0.0, scale * -0.40451, scale * -0.29389);
    rModelPart.CreateNewNode(initial_node_id + 27, scale * 0.47553, scale * -0.15451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 28, scale * -0.43653, scale * 0.18066, scale * -0.16372);
    rModelPart.CreateNewNode(initial_node_id + 29, scale * 0.43653, scale * 0.18066, scale * 0.16372);
    rModelPart.CreateNewNode(initial_node_id + 30, scale * -0.29389, scale * 0.40451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 31, 0.0, scale * 0.47553, scale * 0.15451);
    rModelPart.CreateNewNode(initial_node_id + 32, scale * 0.40451, scale * -0.29389, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 33, scale * 0.23141, scale * 0.37802, scale * 0.23141);
    rModelPart.CreateNewNode(initial_node_id + 34, scale * -0.23141, scale * -0.37802, scale * -0.23141);
    rModelPart.CreateNewNode(initial_node_id + 35, scale * -0.40451, scale * 0.29389, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 36, 0.0, scale * -0.47553, scale * -0.15451);
    rModelPart.CreateNewNode(initial_node_id + 37, scale * 0.29389, scale * -0.40451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 38, scale * -0.43653, scale * -0.18066, scale * -0.16372);
    rModelPart.CreateNewNode(initial_node_id + 39, scale * 0.43653, scale * -0.18066, scale * 0.16372);
    rModelPart.CreateNewNode(initial_node_id + 40, scale * -0.47553, scale * 0.15451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 41, 0.0, scale * 0.40451, scale * 0.29389);
    rModelPart.CreateNewNode(initial_node_id + 42, scale * 0.15451, scale * -0.47553, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 43, scale * -0.23141, scale * 0.37802, scale * 0.23141);
    rModelPart.CreateNewNode(initial_node_id + 44, scale * -0.5, scale * -0, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 45, scale * 0.36663, scale * -0, scale * 0.33997);
    rModelPart.CreateNewNode(initial_node_id + 46, 0.0, scale * -0.5, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 47, scale * -0.47553, scale * -0.15451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 48, scale * -0.43653, scale * 0.18066, scale * 0.16372);
    rModelPart.CreateNewNode(initial_node_id + 49, scale * -0.15451, scale * -0.47553, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 50, scale * -0.40451, scale * -0.29389, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 51, scale * -0.29389, scale * -0.40451, 0.0);
    rModelPart.CreateNewNode(initial_node_id + 52, 0.0, scale * 0.29389, scale * 0.40451);
    rModelPart.CreateNewNode(initial_node_id + 53, scale * 0.16372, scale * 0.18066, scale * 0.43653);
    rModelPart.CreateNewNode(initial_node_id + 54, scale * 0.23141, scale * -0.37802, scale * 0.23141);
    rModelPart.CreateNewNode(initial_node_id + 55, 0.0, scale * -0.47553, scale * 0.15451);
    rModelPart.CreateNewNode(initial_node_id + 56, scale * -0.43653, scale * -0.18066, scale * 0.16372);
    rModelPart.CreateNewNode(initial_node_id + 57, 0.0, scale * 0.15451, scale * 0.47553);
    rModelPart.CreateNewNode(initial_node_id + 58, scale * -0.16372, scale * 0.18066, scale * 0.43653);
    rModelPart.CreateNewNode(initial_node_id + 59, scale * 0.16372, scale * -0.18066, scale * 0.43653);
    rModelPart.CreateNewNode(initial_node_id + 60, 0.0, scale * -0.40451, scale * 0.29389);
    rModelPart.CreateNewNode(initial_node_id + 61, scale * -0.23141, scale * -0.37802, scale * 0.23141);
    rModelPart.CreateNewNode(initial_node_id + 62, scale * -0.33997, scale * -0, scale * 0.36663);
    rModelPart.CreateNewNode(initial_node_id + 63, 0.0, scale * -0, scale * 0.5);
    rModelPart.CreateNewNode(initial_node_id + 64, 0.0, scale * -0.29389, scale * 0.40451);
    rModelPart.CreateNewNode(initial_node_id + 65, 0.0, scale * -0.15451, scale * 0.47553);
    rModelPart.CreateNewNode(initial_node_id + 66, scale * -0.16372, scale * -0.18066, scale * 0.43653);

    // Modify center
    block_for_each(rModelPart.Nodes(), [&rCenter](Node& rNode) {
        rNode.X()  += rCenter[0];
        rNode.X0() += rCenter[0];
        rNode.Y()  += rCenter[1];
        rNode.Y0() += rCenter[1];
        rNode.Z()  += rCenter[2];
        rNode.Z0() += rCenter[2];
    });

    // Get the Initial Id
    const std::size_t initial_condition_id = block_for_each<MaxReduction<std::size_t>>(rModelPart.GetRootModelPart().Conditions(), [](Condition& rCond) {
        return rCond.Id();
    });

    // Now we create the conditions
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 1,   {{initial_node_id + 53,  initial_node_id + 29, initial_node_id + 33}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 2,   {{initial_node_id + 8,   initial_node_id + 28, initial_node_id + 13}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 3,   {{initial_node_id + 66,  initial_node_id + 62, initial_node_id + 56}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 4,   {{initial_node_id + 48,  initial_node_id + 56, initial_node_id + 62}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 5,   {{initial_node_id + 62,  initial_node_id + 58, initial_node_id + 48}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 6,   {{initial_node_id + 62,  initial_node_id + 66, initial_node_id + 58}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 7,   {{initial_node_id + 22,  initial_node_id + 38, initial_node_id + 28}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 8,   {{initial_node_id + 39,  initial_node_id + 45, initial_node_id + 59}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 9,   {{initial_node_id + 45,  initial_node_id + 29, initial_node_id + 53}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 10,  {{initial_node_id + 45,  initial_node_id + 39, initial_node_id + 29}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 11,  {{initial_node_id + 5,   initial_node_id + 9,  initial_node_id + 1}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 12,  {{initial_node_id + 39,  initial_node_id + 59, initial_node_id + 54}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 13,  {{initial_node_id + 9,   initial_node_id + 24, initial_node_id + 26}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 14,  {{initial_node_id + 39,  initial_node_id + 54, initial_node_id + 37}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 15,  {{initial_node_id + 48,  initial_node_id + 43, initial_node_id + 30}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 16,  {{initial_node_id + 6,   initial_node_id + 11, initial_node_id + 1}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 17,  {{initial_node_id + 8,   initial_node_id + 13, initial_node_id + 7}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 18,  {{initial_node_id + 38,  initial_node_id + 22, initial_node_id + 14}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 19,  {{initial_node_id + 8,   initial_node_id + 14, initial_node_id + 22}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 20,  {{initial_node_id + 26,  initial_node_id + 34, initial_node_id + 14}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 21,  {{initial_node_id + 61,  initial_node_id + 66, initial_node_id + 56}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 22,  {{initial_node_id + 37,  initial_node_id + 24, initial_node_id + 19}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 23,  {{initial_node_id + 7,   initial_node_id + 6,  initial_node_id + 1}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 24,  {{initial_node_id + 41,  initial_node_id + 53, initial_node_id + 33}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 25,  {{initial_node_id + 50,  initial_node_id + 38, initial_node_id + 51}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 26,  {{initial_node_id + 32,  initial_node_id + 39, initial_node_id + 37}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 27,  {{initial_node_id + 3,   initial_node_id + 8,  initial_node_id + 7}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 28,  {{initial_node_id + 15,  initial_node_id + 9,  initial_node_id + 26}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 29,  {{initial_node_id + 58,  initial_node_id + 52, initial_node_id + 41}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 30,  {{initial_node_id + 52,  initial_node_id + 53, initial_node_id + 41}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 31,  {{initial_node_id + 36,  initial_node_id + 34, initial_node_id + 26}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 32,  {{initial_node_id + 29,  initial_node_id + 17, initial_node_id + 16}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 33,  {{initial_node_id + 58,  initial_node_id + 66, initial_node_id + 63}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 34,  {{initial_node_id + 63,  initial_node_id + 66, initial_node_id + 65}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 35,  {{initial_node_id + 4,   initial_node_id + 8,  initial_node_id + 2}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 36,  {{initial_node_id + 56,  initial_node_id + 50, initial_node_id + 51}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 37,  {{initial_node_id + 23,  initial_node_id + 11, initial_node_id + 20}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 38,  {{initial_node_id + 47,  initial_node_id + 56, initial_node_id + 44}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 39,  {{initial_node_id + 29,  initial_node_id + 39, initial_node_id + 23}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 40,  {{initial_node_id + 23,  initial_node_id + 39, initial_node_id + 27}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 41,  {{initial_node_id + 42,  initial_node_id + 36, initial_node_id + 24}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 42,  {{initial_node_id + 65,  initial_node_id + 59, initial_node_id + 63}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 43,  {{initial_node_id + 48,  initial_node_id + 58, initial_node_id + 43}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 44,  {{initial_node_id + 57,  initial_node_id + 58, initial_node_id + 63}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 45,  {{initial_node_id + 50,  initial_node_id + 47, initial_node_id + 38}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 46,  {{initial_node_id + 35,  initial_node_id + 48, initial_node_id + 30}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 47,  {{initial_node_id + 40,  initial_node_id + 28, initial_node_id + 44}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 48,  {{initial_node_id + 19,  initial_node_id + 32, initial_node_id + 37}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 49,  {{initial_node_id + 2,   initial_node_id + 1,  initial_node_id + 4}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 50,  {{initial_node_id + 5,   initial_node_id + 1,  initial_node_id + 11}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 51,  {{initial_node_id + 10,  initial_node_id + 14, initial_node_id + 4}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 52,  {{initial_node_id + 14,  initial_node_id + 15, initial_node_id + 26}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 53,  {{initial_node_id + 46,  initial_node_id + 55, initial_node_id + 49}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 54,  {{initial_node_id + 3,   initial_node_id + 1,  initial_node_id + 2}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 55,  {{initial_node_id + 42,  initial_node_id + 24, initial_node_id + 37}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 56,  {{initial_node_id + 11,  initial_node_id + 6,  initial_node_id + 16}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 57,  {{initial_node_id + 50,  initial_node_id + 56, initial_node_id + 47}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 58,  {{initial_node_id + 18,  initial_node_id + 12, initial_node_id + 21}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 59,  {{initial_node_id + 27,  initial_node_id + 19, initial_node_id + 23}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 60,  {{initial_node_id + 15,  initial_node_id + 14, initial_node_id + 10}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 61,  {{initial_node_id + 18,  initial_node_id + 21, initial_node_id + 31}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 62,  {{initial_node_id + 32,  initial_node_id + 19, initial_node_id + 27}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 63,  {{initial_node_id + 12,  initial_node_id + 18, initial_node_id + 6}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 64,  {{initial_node_id + 21,  initial_node_id + 12, initial_node_id + 25}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 65,  {{initial_node_id + 34,  initial_node_id + 38, initial_node_id + 14}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 66,  {{initial_node_id + 46,  initial_node_id + 36, initial_node_id + 42}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 67,  {{initial_node_id + 35,  initial_node_id + 28, initial_node_id + 40}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 68,  {{initial_node_id + 28,  initial_node_id + 35, initial_node_id + 30}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 69,  {{initial_node_id + 64,  initial_node_id + 65, initial_node_id + 66}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 70,  {{initial_node_id + 5,   initial_node_id + 19, initial_node_id + 9}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 71,  {{initial_node_id + 20,  initial_node_id + 29, initial_node_id + 23}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 72,  {{initial_node_id + 1,   initial_node_id + 3,  initial_node_id + 7}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 73,  {{initial_node_id + 63,  initial_node_id + 53, initial_node_id + 57}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 74,  {{initial_node_id + 1,   initial_node_id + 9,  initial_node_id + 4}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 75,  {{initial_node_id + 19,  initial_node_id + 5,  initial_node_id + 11}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 76,  {{initial_node_id + 17,  initial_node_id + 29, initial_node_id + 20}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 77,  {{initial_node_id + 4,   initial_node_id + 9,  initial_node_id + 10}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 78,  {{initial_node_id + 46,  initial_node_id + 42, initial_node_id + 55}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 79,  {{initial_node_id + 35,  initial_node_id + 40, initial_node_id + 48}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 80,  {{initial_node_id + 52,  initial_node_id + 58, initial_node_id + 57}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 81,  {{initial_node_id + 31,  initial_node_id + 41, initial_node_id + 33}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 82,  {{initial_node_id + 66,  initial_node_id + 61, initial_node_id + 60}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 83,  {{initial_node_id + 59,  initial_node_id + 64, initial_node_id + 60}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 84,  {{initial_node_id + 52,  initial_node_id + 57, initial_node_id + 53}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 85,  {{initial_node_id + 31,  initial_node_id + 43, initial_node_id + 41}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 86,  {{initial_node_id + 36,  initial_node_id + 26, initial_node_id + 24}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 87,  {{initial_node_id + 64,  initial_node_id + 59, initial_node_id + 65}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 88,  {{initial_node_id + 32,  initial_node_id + 27, initial_node_id + 39}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 89,  {{initial_node_id + 12,  initial_node_id + 13, initial_node_id + 25}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 90,  {{initial_node_id + 42,  initial_node_id + 54, initial_node_id + 55}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 91,  {{initial_node_id + 31,  initial_node_id + 33, initial_node_id + 18}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 92,  {{initial_node_id + 12,  initial_node_id + 7,  initial_node_id + 13}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 93,  {{initial_node_id + 22,  initial_node_id + 28, initial_node_id + 8}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 94,  {{initial_node_id + 64,  initial_node_id + 66, initial_node_id + 60}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 95,  {{initial_node_id + 3,   initial_node_id + 2,  initial_node_id + 8}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 96,  {{initial_node_id + 49,  initial_node_id + 34, initial_node_id + 36}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 97,  {{initial_node_id + 18,  initial_node_id + 33, initial_node_id + 16}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 98,  {{initial_node_id + 42,  initial_node_id + 37, initial_node_id + 54}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 99,  {{initial_node_id + 55,  initial_node_id + 54, initial_node_id + 60}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 100, {{initial_node_id + 45,  initial_node_id + 53, initial_node_id + 59}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 101, {{initial_node_id + 49,  initial_node_id + 61, initial_node_id + 51}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 102, {{initial_node_id + 24,  initial_node_id + 9,  initial_node_id + 19}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 103, {{initial_node_id + 25,  initial_node_id + 13, initial_node_id + 30}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 104, {{initial_node_id + 49,  initial_node_id + 51, initial_node_id + 34}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 105, {{initial_node_id + 12,  initial_node_id + 6,  initial_node_id + 7}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 106, {{initial_node_id + 30,  initial_node_id + 13, initial_node_id + 28}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 107, {{initial_node_id + 59,  initial_node_id + 53, initial_node_id + 63}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 108, {{initial_node_id + 44,  initial_node_id + 48, initial_node_id + 40}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 109, {{initial_node_id + 18,  initial_node_id + 16, initial_node_id + 6}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 110, {{initial_node_id + 25,  initial_node_id + 31, initial_node_id + 21}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 111, {{initial_node_id + 19,  initial_node_id + 11, initial_node_id + 23}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 112, {{initial_node_id + 28,  initial_node_id + 38, initial_node_id + 44}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 113, {{initial_node_id + 17,  initial_node_id + 11, initial_node_id + 16}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 114, {{initial_node_id + 15,  initial_node_id + 10, initial_node_id + 9}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 115, {{initial_node_id + 14,  initial_node_id + 8,  initial_node_id + 4}},  p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 116, {{initial_node_id + 16,  initial_node_id + 33, initial_node_id + 29}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 117, {{initial_node_id + 60,  initial_node_id + 54, initial_node_id + 59}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 118, {{initial_node_id + 51,  initial_node_id + 61, initial_node_id + 56}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 119, {{initial_node_id + 41,  initial_node_id + 43, initial_node_id + 58}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 120, {{initial_node_id + 36,  initial_node_id + 46, initial_node_id + 49}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 121, {{initial_node_id + 17,  initial_node_id + 20, initial_node_id + 11}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 122, {{initial_node_id + 25,  initial_node_id + 30, initial_node_id + 43}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 123, {{initial_node_id + 49,  initial_node_id + 55, initial_node_id + 61}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 124, {{initial_node_id + 56,  initial_node_id + 48, initial_node_id + 44}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 125, {{initial_node_id + 55,  initial_node_id + 60, initial_node_id + 61}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 126, {{initial_node_id + 31,  initial_node_id + 25, initial_node_id + 43}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 127, {{initial_node_id + 44,  initial_node_id + 38, initial_node_id + 47}}, p_prop);
    rModelPart.CreateNewCondition(rConditionName, initial_condition_id + 128, {{initial_node_id + 51,  initial_node_id + 38, initial_node_id + 34}}, p_prop);
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& CreateCubeSkinModelPart(
    Model& rCurrentModel,
    const double HalfX,
    const double HalfY,
    const double HalfZ,
    const DataCommunicator& rDataCommunicator
    )
{
    // Generate the cube skin
    ModelPart& r_skin_part = rCurrentModel.CreateModelPart("Skin");

    // Distributed related variables
    if (rDataCommunicator.IsDistributed()) {
        r_skin_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    }
    const int rank =  rDataCommunicator.Rank();
    const int world_size = rDataCommunicator.Size();

    // Create properties
    auto p_properties = r_skin_part.CreateNewProperties(1, 0);

    // If only one partition
    if (world_size == 1) {
        // Create nodes
        r_skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
        r_skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
        r_skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
        r_skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);
        r_skin_part.CreateNewNode(5, -HalfX, -HalfY,  HalfZ);
        r_skin_part.CreateNewNode(6,  HalfX, -HalfY,  HalfZ);
        r_skin_part.CreateNewNode(7,  HalfX,  HalfY,  HalfZ);
        r_skin_part.CreateNewNode(8, -HalfX,  HalfY,  HalfZ);

        // Set the partition index
        if (rDataCommunicator.IsDistributed()) {
            for (auto& r_node : r_skin_part.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            }
        }

        // Create elements
        r_skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
        r_skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);
    } else {
        if (rank == 0) {
            // Create nodes
            auto p_node1 = r_skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
            auto p_node2 = r_skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
            auto p_node3 = r_skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
            auto p_node4 = r_skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);

            // Set partitions
            p_node1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node2->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node3->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node4->FastGetSolutionStepValue(PARTITION_INDEX) = 0;

            // Create elements
            r_skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
        } else if (rank == 1) {
            // Create nodes
            auto p_node1 = r_skin_part.CreateNewNode(1, -HalfX, -HalfY, -HalfZ);
            auto p_node2 = r_skin_part.CreateNewNode(2,  HalfX, -HalfY, -HalfZ);
            auto p_node3 = r_skin_part.CreateNewNode(3,  HalfX,  HalfY, -HalfZ);
            auto p_node4 = r_skin_part.CreateNewNode(4, -HalfX,  HalfY, -HalfZ);
            auto p_node5 = r_skin_part.CreateNewNode(5, -HalfX, -HalfY,  HalfZ);
            auto p_node6 = r_skin_part.CreateNewNode(6,  HalfX, -HalfY,  HalfZ);
            auto p_node7 = r_skin_part.CreateNewNode(7,  HalfX,  HalfY,  HalfZ);
            auto p_node8 = r_skin_part.CreateNewNode(8, -HalfX,  HalfY,  HalfZ);

            // Set partitions
            p_node1->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node2->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node3->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node4->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node5->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            p_node6->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            p_node7->FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            p_node8->FastGetSolutionStepValue(PARTITION_INDEX) = 1;

            // Create elements
            r_skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
            r_skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);
        }
    }

    return r_skin_part;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& CreateCubeModelPart(
    Model& rCurrentModel,
    const DataCommunicator& rDataCommunicator
    )
{
    // Generate the cube skin
    ModelPart& r_model_part = rCurrentModel.CreateModelPart("Cube");

    // Distributed related variables
    if (rDataCommunicator.IsDistributed()) {
        r_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    }
    const int rank =  rDataCommunicator.Rank();
    const int world_size = rDataCommunicator.Size();

    // Create properties
    auto p_properties = r_model_part.CreateNewProperties(1, 0);

    if (world_size == 1) {
        // Create nodes
        r_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
        r_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
        r_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
        r_model_part.CreateNewNode(4 , 0.5 , 1.0 , 1.0);
        r_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
        r_model_part.CreateNewNode(6 , 0.5 , 1.0 , 0.0);
        r_model_part.CreateNewNode(7 , 0.5 , 0.0 , 1.0);
        r_model_part.CreateNewNode(8 , 0.5 , 0.0 , 0.0);
        r_model_part.CreateNewNode(9 , 1.0 , 1.0 , 1.0);
        r_model_part.CreateNewNode(10 , 1.0 , 1.0 , 0.0);
        r_model_part.CreateNewNode(11 , 1.0 , 0.0 , 1.0);
        r_model_part.CreateNewNode(12 , 1.0 , 0.0 , 0.0);

        // Set the partition index
        if (rDataCommunicator.IsDistributed()) {
            for (auto& r_node : r_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            }
        }

        // Create elements
        r_model_part.CreateNewElement("Element3D8N", 1, {{5,8,6,2,3,7,4,1}}, p_properties);
        r_model_part.CreateNewElement("Element3D8N", 2, {{8,12,10,6,7,11,9,4}}, p_properties);
    } else { // Assuming always two partitions
        if (rank == 0) {
            // Create nodes
            r_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            r_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            r_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            r_model_part.CreateNewNode(4 , 0.5 , 1.0 , 1.0);
            r_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            r_model_part.CreateNewNode(6 , 0.5 , 1.0 , 0.0);
            r_model_part.CreateNewNode(7 , 0.5 , 0.0 , 1.0);
            r_model_part.CreateNewNode(8 , 0.5 , 0.0 , 0.0);

            // Set partitions
            for (auto& r_node : r_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            }

            // Create elements
            r_model_part.CreateNewElement("Element3D8N", 1, {{5,8,6,2,3,7,4,1}}, p_properties);
        } else if (rank == 1) {
            // Create nodes
            auto p_node4  = r_model_part.CreateNewNode(4 , 0.5 , 1.0 , 1.0);
            auto p_node6  = r_model_part.CreateNewNode(6 , 0.5 , 1.0 , 0.0);
            auto p_node7  = r_model_part.CreateNewNode(7 , 0.5 , 0.0 , 1.0);
            auto p_node8  = r_model_part.CreateNewNode(8 , 0.5 , 0.0 , 0.0);
            auto p_node9  = r_model_part.CreateNewNode(9 , 1.0 , 1.0 , 1.0);
            auto p_node10 = r_model_part.CreateNewNode(10 , 1.0 , 1.0 , 0.0);
            auto p_node11 = r_model_part.CreateNewNode(11 , 1.0 , 0.0 , 1.0);
            auto p_node12 = r_model_part.CreateNewNode(12 , 1.0 , 0.0 , 0.0);

            // Set partitions
            for (auto& r_node : r_model_part.Nodes()) {
                r_node.FastGetSolutionStepValue(PARTITION_INDEX) = 1;
            }
            p_node4->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node6->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node7->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
            p_node8->FastGetSolutionStepValue(PARTITION_INDEX) = 0;

            // Create elements
            r_model_part.CreateNewElement("Element3D8N", 2, {{8,12,10,6,7,11,9,4}}, p_properties);
        }
    }

    return r_model_part;
}

} // namespace ConstraintUtilities
} // namespace Kratos