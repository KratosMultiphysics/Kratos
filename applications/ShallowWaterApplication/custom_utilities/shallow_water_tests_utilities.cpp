//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes

// External includes

// Project includes
#include "shallow_water_tests_utilities.h"

namespace Kratos {

void ShallowWaterTestsUtilities::AssembleRHS(
    Vector& rRHS_element,
    const Vector& rRHS_condition,
    const std::vector<IndexType>& rIds)
{
    IndexType n_dofs = rRHS_condition.size() / rIds.size();
    for (IndexType i = 0; i < rRHS_condition.size(); ++i)
    {
        IndexType i_dof = i % n_dofs;
        IndexType local_id = i / n_dofs;
        IndexType global_id = rIds[local_id] - 1;
        IndexType elem_pos = i_dof + global_id * n_dofs;
        rRHS_element[elem_pos] += rRHS_condition[i];
    }
}

void ShallowWaterTestsUtilities::AddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(MOMENTUM);
    rModelPart.AddNodalSolutionStepVariable(HEIGHT);
    rModelPart.AddNodalSolutionStepVariable(TOPOGRAPHY);
    rModelPart.AddNodalSolutionStepVariable(RAIN);
    rModelPart.AddNodalSolutionStepVariable(MANNING);
    rModelPart.AddNodalSolutionStepVariable(WIND);
    rModelPart.AddNodalSolutionStepVariable(ATMOSPHERIC_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(VERTICAL_VELOCITY);
}

void TestCreateGeometry(
    ModelPart& rModelPart,
    std::string ElementName,
    std::string ConditionName)
{
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer property = rModelPart.CreateNewProperties(0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    rModelPart.CreateNewElement(ElementName, 1, elem_nodes, property);
    std::vector<IndexType> cond_1_nodes {1, 2};
    std::vector<IndexType> cond_2_nodes {2, 3};
    std::vector<IndexType> cond_3_nodes {3, 1};
    rModelPart.CreateNewCondition(ConditionName, 1, cond_1_nodes, property);
    rModelPart.CreateNewCondition(ConditionName, 2, cond_2_nodes, property);
    rModelPart.CreateNewCondition(ConditionName, 3, cond_3_nodes, property);
}

} // namespace Kratos
