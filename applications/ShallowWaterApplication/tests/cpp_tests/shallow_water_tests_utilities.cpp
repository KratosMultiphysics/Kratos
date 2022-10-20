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
#include "shallow_water_application_variables.h"

namespace Kratos {

namespace Testing {

void ShallowWaterTestsUtilities::AssembleRHS(
    Vector& rRHS_element,
    const Vector& rRHS_condition,
    const std::vector<std::size_t>& rIds)
{
    std::size_t n_dofs = rRHS_condition.size() / rIds.size();
    for (std::size_t i = 0; i < rRHS_condition.size(); ++i)
    {
        std::size_t i_dof = i % n_dofs;
        std::size_t local_id = i / n_dofs;
        std::size_t global_id = rIds[local_id] - 1;
        std::size_t elem_pos = i_dof + global_id * n_dofs;
        rRHS_element[elem_pos] += rRHS_condition[i];
    }
}

void ShallowWaterTestsUtilities::AddVariables(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(MOMENTUM);
    rModelPart.AddNodalSolutionStepVariable(HEIGHT);
    rModelPart.AddNodalSolutionStepVariable(FREE_SURFACE_ELEVATION);
    rModelPart.AddNodalSolutionStepVariable(TOPOGRAPHY);
    rModelPart.AddNodalSolutionStepVariable(RAIN);
    rModelPart.AddNodalSolutionStepVariable(MANNING);
    rModelPart.AddNodalSolutionStepVariable(WIND);
    rModelPart.AddNodalSolutionStepVariable(ATMOSPHERIC_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(VERTICAL_VELOCITY);
}

void ShallowWaterTestsUtilities::CreateGeometry(
    ModelPart& rModelPart,
    const std::string& rElementName,
    const std::string& rConditionName)
{
    Properties::Pointer property = rModelPart.CreateNewProperties(0);

    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

    rModelPart.CreateNewElement(rElementName, 1, {{1, 2, 3}}, property);

    rModelPart.CreateNewCondition(rConditionName, 1, {{1, 2}}, property);
    rModelPart.CreateNewCondition(rConditionName, 2, {{2, 3}}, property);
    rModelPart.CreateNewCondition(rConditionName, 3, {{3, 1}}, property);
}

void ShallowWaterTestsUtilities::CalculateAndAssembleRHS(
    ModelPart& rModelPart,
    Vector& rRHS)
{
    Matrix lhs = ZeroMatrix(9,9);
    Vector rhs_cond = ZeroVector(6);
    auto& process_info = rModelPart.GetProcessInfo();

    auto& element = rModelPart.GetElement(1);
    element.CalculateLocalSystem(lhs, rRHS, process_info);

    auto& condition_1 = rModelPart.GetCondition(1);
    condition_1.CalculateRightHandSide(rhs_cond, process_info);
    AssembleRHS(rRHS, rhs_cond, {1, 2});

    auto& condition_2 = rModelPart.GetCondition(2);
    condition_2.CalculateRightHandSide(rhs_cond, process_info);
    AssembleRHS(rRHS, rhs_cond, {2, 3});

    auto& condition_3 = rModelPart.GetCondition(3);
    condition_3.CalculateRightHandSide(rhs_cond, process_info);
    AssembleRHS(rRHS, rhs_cond, {3, 1});
}

} // namespace Testing

} // namespace Kratos
