//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Vicente Mataix Ferrandiz
//

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"

#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "processes/skin_detection_process.h"
#include "includes/element.h"

namespace Kratos
{
namespace Testing
{

typedef Node<3> NodeType;

KRATOS_TEST_CASE_IN_SUITE(SkinDetectionProcess, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("test_model_part",2);

    model_part.AddNodalSolutionStepVariable(TEMPERATURE);;

    NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = model_part.CreateNewNode(4, 0.0, 1.0, 0.0);

    Geometry<NodeType>::PointsArrayType array_nodes1, array_nodes2;
    array_nodes1.push_back(p_node_1);
    array_nodes1.push_back(p_node_2);
    array_nodes1.push_back(p_node_3);

    array_nodes2.push_back(p_node_1);
    array_nodes2.push_back(p_node_3);
    array_nodes2.push_back(p_node_4);

    Properties::Pointer p_elem_prop = model_part.pGetProperties(0);
    Element::Pointer p_elem_1 = model_part.CreateNewElement("Element2D3N", 1,PointerVector<NodeType>{array_nodes1}, p_elem_prop);
    Element::Pointer p_elem_2 = model_part.CreateNewElement("Element2D3N", 2,PointerVector<NodeType>{array_nodes2}, p_elem_prop);

    Parameters default_parameters = Parameters(R"(
        {
            "name_auxiliar_model_part"              : "SkinModelPart",
            "name_auxiliar_condition"               : "Condition",
            "list_model_parts_to_assign_conditions" : [],
            "echo_level"                            : 0
        })");

    SkinDetectionProcess<2> skin_process = SkinDetectionProcess<2>(
                                           model_part, default_parameters);

    // We generate in several iterations to see if it crashes
    for (int i = 0; i < 2; i++) {
        skin_process.Execute();
        KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("SkinModelPart").NumberOfConditions(), 4);
    }

    // Now we remove one element
    p_elem_2->Set(TO_ERASE);
    model_part.RemoveElementsFromAllLevels(TO_ERASE);

    // We execute again
    skin_process.Execute();
    KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("SkinModelPart").NumberOfConditions(), 3);
}
} // namespace Testing
} // namespace Kratos
