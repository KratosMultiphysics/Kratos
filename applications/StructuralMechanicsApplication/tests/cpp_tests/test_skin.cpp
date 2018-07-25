// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
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

KRATOS_TEST_CASE_IN_SUITE(SkinDetectionProcess, KratosStructuralMechanicsFastSuite)
{
    ModelPart model_part("test_model_part");

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

    Properties::Pointer material_properties;
    material_properties->SetValue(YOUNG_MODULUS, 210e9);
    material_properties->SetValue(POISSON_RATIO, 0.22);

    Element::Pointer pElem1 = Element().Clone(1,array_nodes1);
    Element::Pointer pElem2 = Element().Clone(1,array_nodes2);

    model_part.AddNode(p_node_1);
    model_part.AddNode(p_node_2);
    model_part.AddNode(p_node_3);
    model_part.AddNode(p_node_4);

    model_part.AddElement(pElem1);
    model_part.AddElement(pElem2);

    Parameters default_parameters = Parameters(R"(
        {
            "name_auxiliar_model_part"              : "SkinModelPart",
            "name_auxiliar_condition"               : "Condition",
            "list_model_parts_to_assign_conditions" : [],
            "echo_level"                            : 0
        })");

    SkinDetectionProcess<2> skin_process = SkinDetectionProcess<2>(
                                           model_part, default_parameters);

    for (int i = 0; i < 2; i++)
    {
        
    }
}

} // namespace Testing
} // namespace Kratos