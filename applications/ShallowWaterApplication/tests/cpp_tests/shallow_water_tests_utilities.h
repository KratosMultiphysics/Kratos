//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "shallow_water_application_variables.h"

void AssembleRHS_Local(Vector& rRHS_element, const Vector& rRHS_condition, const std::vector<IndexType>& rIds)
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

void AddVariables(ModelPart& rModelPart)
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
    rModelPart.AddNodalSolutionStepVariable(MESH_ACCELERATION);
}

void CreateGeometry(ModelPart$ rModelPart, std::string ElementName, std::string ConditionName)
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

void SetNodalValues(
    ModelPart& rModelPart,
    const double& rManning,
    const double& rHeight,
    const array_1d<double,3>& rVelocity,
    const array_1d<double,3>& rTopographySlope,
    const array_1d<double,3>& rHeightGradient)
{
    for (auto& r_node : rModelPart.Nodes())
    {
        const array_1d<double,3> coords = r_node.Coordinates();
        const auto height = rHeight + inner_prod(coords, rHeightGradient);
        const auto velocity = rMomentum / height;
        const auto topography = inner_prod(coords, rTopographySlope);

        r_node.FastGetSolutionStepValue(HEIGHT) = height;
        r_node.FastGetSolutionStepValue(VELOCITY) = velocity;
        r_node.FastGetSolutionStepValue(MOMENTUM) = rMomentum;
        r_node.FastGetSolutionStepValue(MANNING) = rManning;
        r_node.FastGetSolutionStepValue(TOPOGRAPHY) = topography;
    }
}

void AssembleRHS(ModelPart& rModelPart, Vector& rRHS)
{
    // Compute RHS and LHS for the element
    rRhs = ZeroVector(9);
    Matrix lhs = ZeroMatrix(9,9);
    const ProcessInfo& r_process_info = model_part.GetProcessInfo();
    rModelPart.GetElement(1).CalculateLocalSystem(lhs, rRhs, r_process_info);

    // Add the conditions contributions
    Vector rhs_cond = ZeroVector(6);
    rModelPart.GetCondition(1).CalculateRightHandSide(rhs_cond, r_process_info);
    AssembleRHS_Local(rRhs, rhs_cond, {1, 2});

    rModelPart.GetCondition(2).CalculateRightHandSide(rhs_cond, r_process_info);
    AssembleRHS_Local(rRhs, rhs_cond, {2, 3});

    rModelPart.GetCondition(3).CalculateRightHandSide(rhs_cond, r_process_info);
    AssembleRHS_Local(rRhs, rhs_cond, {3, 1});
}
