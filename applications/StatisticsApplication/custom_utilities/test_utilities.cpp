//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

// Application incl

// Include base h
#include "test_utilities.h"

namespace Kratos
{
namespace StatisticsApplicationTestUtilities
{
using NodeType = ModelPart::NodeType;
using NodesContainerType = ModelPart::NodesContainerType;
using ConditionType = ModelPart::ConditionType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;
using ElementType = ModelPart::ElementType;
using ElementsContainerType = ModelPart::ElementsContainerType;

void AddNodalSolutionStepVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
}

void CreateModelPart(ModelPart& rModelPart)
{
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(4, 2.0, 1.0, 0.0);
    rModelPart.CreateNewNode(5, 1.5, 1.0, 0.0);
    rModelPart.CreateNewNode(6, 0.5, 1.0, 0.0);
    rModelPart.CreateNewNode(7, 0.0, 1.0, 0.0);

    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

    std::vector<ModelPart::IndexType> elem_1_nodes{1, 6, 7};
    std::vector<ModelPart::IndexType> elem_2_nodes{1, 2, 6};
    std::vector<ModelPart::IndexType> elem_3_nodes{6, 2, 5};
    std::vector<ModelPart::IndexType> elem_4_nodes{2, 3, 5};
    std::vector<ModelPart::IndexType> elem_5_nodes{5, 3, 4};

    rModelPart.CreateNewElement("Element2D3N",
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_1_nodes, p_elem_prop);
    rModelPart.CreateNewElement("Element2D3N",
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_2_nodes, p_elem_prop);
    rModelPart.CreateNewElement("Element2D3N",
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_3_nodes, p_elem_prop);
    rModelPart.CreateNewElement("Element2D3N",
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_4_nodes, p_elem_prop);
    rModelPart.CreateNewElement("Element2D3N",
                                rModelPart.GetRootModelPart().NumberOfElements() + 1,
                                elem_5_nodes, p_elem_prop);

    std::vector<ModelPart::IndexType> cond_1_nodes{1, 2};
    std::vector<ModelPart::IndexType> cond_2_nodes{2, 3};
    std::vector<ModelPart::IndexType> cond_3_nodes{3, 4};
    std::vector<ModelPart::IndexType> cond_4_nodes{4, 5};
    std::vector<ModelPart::IndexType> cond_5_nodes{5, 6};
    std::vector<ModelPart::IndexType> cond_6_nodes{6, 7};
    std::vector<ModelPart::IndexType> cond_7_nodes{7, 1};

    rModelPart.CreateNewCondition(
        "Condition2D2N", rModelPart.GetRootModelPart().NumberOfConditions() + 1,
        cond_1_nodes, p_elem_prop);
    rModelPart.CreateNewCondition(
        "Condition2D2N", rModelPart.GetRootModelPart().NumberOfConditions() + 1,
        cond_2_nodes, p_elem_prop);
    rModelPart.CreateNewCondition(
        "Condition2D2N", rModelPart.GetRootModelPart().NumberOfConditions() + 1,
        cond_3_nodes, p_elem_prop);
    rModelPart.CreateNewCondition(
        "Condition2D2N", rModelPart.GetRootModelPart().NumberOfConditions() + 1,
        cond_4_nodes, p_elem_prop);
    rModelPart.CreateNewCondition(
        "Condition2D2N", rModelPart.GetRootModelPart().NumberOfConditions() + 1,
        cond_5_nodes, p_elem_prop);
    rModelPart.CreateNewCondition(
        "Condition2D2N", rModelPart.GetRootModelPart().NumberOfConditions() + 1,
        cond_6_nodes, p_elem_prop);
    rModelPart.CreateNewCondition(
        "Condition2D2N", rModelPart.GetRootModelPart().NumberOfConditions() + 1,
        cond_7_nodes, p_elem_prop);
}

} // namespace StatisticsApplicationTestUtilities
} // namespace Kratos