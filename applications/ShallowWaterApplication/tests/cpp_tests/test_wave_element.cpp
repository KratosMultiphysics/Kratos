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
#include "testing/testing.h"
#include "containers/model.h"
#include "shallow_water_application_variables.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType IndexType;

void AssembleRHS_condition(Vector& rRHS_element, const Vector& rRHS_condition, const std::vector<IndexType>& rIds)
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

void WaveElementSteadyStateTest(
    const double& rManning,
    const double& rHeight,
    const array_1d<double,3>& rVelocity,
    const array_1d<double,3>& rTopographySlope,
    const array_1d<double,3>& rHeightGradient = ZeroVector(3),
    const double& rTolerance = 1e-12)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);

    // Variables
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(MANNING);
    model_part.AddNodalSolutionStepVariable(TOPOGRAPHY);

    // Process info creation
    const double gravity = 9.81;
    const double stab_factor = 0.01;
    ProcessInfo& r_process_info = model_part.GetProcessInfo();
    r_process_info.SetValue(GRAVITY_Z, gravity);
    r_process_info.SetValue(STABILIZATION_FACTOR, stab_factor);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer property = model_part.CreateNewProperties(0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    Element::Pointer element = model_part.CreateNewElement("WaveElement2D3N", 1, elem_nodes, property);
    std::vector<IndexType> cond_1_nodes {1, 2};
    std::vector<IndexType> cond_2_nodes {2, 3};
    std::vector<IndexType> cond_3_nodes {3, 1};
    Condition::Pointer cond_1 = model_part.CreateNewCondition("WaveCondition2D2N", 1, cond_1_nodes, property);
    Condition::Pointer cond_2 = model_part.CreateNewCondition("WaveCondition2D2N", 2, cond_2_nodes, property);
    Condition::Pointer cond_3 = model_part.CreateNewCondition("WaveCondition2D2N", 3, cond_3_nodes, property);

    // Set the nodal values
    auto r_geom = element->GetGeometry();
    for (auto& r_node : r_geom)
    {
        const array_1d<double,3> coords = r_node.Coordinates();
        const auto height = rHeight + inner_prod(coords, rHeightGradient);
        const auto topography = inner_prod(coords, rTopographySlope);

        r_node.FastGetSolutionStepValue(VELOCITY) = rVelocity;
        r_node.FastGetSolutionStepValue(HEIGHT) = height;
        r_node.FastGetSolutionStepValue(MANNING) = rManning;
        r_node.FastGetSolutionStepValue(TOPOGRAPHY) = topography;
    }

    // Compute RHS and LHS
    Vector rhs = ZeroVector(9);
    Matrix lhs = ZeroMatrix(9,9);
    element->CalculateLocalSystem(lhs, rhs, r_process_info);

    Vector rhs_cond = ZeroVector(6);
    cond_1->CalculateRightHandSide(rhs_cond, r_process_info);
    AssembleRHS_condition(rhs, rhs_cond, cond_1_nodes);

    cond_2->CalculateRightHandSide(rhs_cond, r_process_info);
    AssembleRHS_condition(rhs, rhs_cond, cond_2_nodes);

    cond_3->CalculateRightHandSide(rhs_cond, r_process_info);
    AssembleRHS_condition(rhs, rhs_cond, cond_3_nodes);

    // Check the RHS values. Since it is a steady solution the RHS must be zero
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(rhs, ZeroVector(9), rTolerance);
}

/**
 * @brief Check the WaveElement2D3N element with still free surface
 */
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3NSteadyStillSurface, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    const array_1d<double,3> velocity = ZeroVector(3);
    const array_1d<double,3> slope = ZeroVector(3);

    WaveElementSteadyStateTest(manning, height, velocity, slope);
}

/**
 * @brief Check the WaveElement2D3N element still free surface and bottom topography x-gradient
 */
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3NSteadyGradient, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    array_1d<double,3> velocity = ZeroVector(3);
    array_1d<double,3> slope = ZeroVector(3);
    slope[0] = 0.05;
    array_1d<double,3> height_gradient = -slope;

    WaveElementSteadyStateTest(manning, height, velocity, slope, height_gradient);
}

/**
 * @brief Check the WaveElement2D3N element still free surface and bottom topography skew gradient
 */
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3NSteadySkewGradient, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    array_1d<double,3> velocity = ZeroVector(3);
    array_1d<double,3> slope = ZeroVector(3);
    slope[0] = 0.05;
    slope[1] = 0.05;
    array_1d<double,3> height_gradient = -slope;

    WaveElementSteadyStateTest(manning, height, velocity, slope, height_gradient);
}

/**
 * @brief Check the WaveElement2D3N element still free surface and bottom topography x-gradient
 */
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3NVelocityAndGradient, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    array_1d<double,3> velocity = ZeroVector(3);
    velocity[0] = 0.6;
    velocity[1] = 1.0;
    array_1d<double,3> slope = ZeroVector(3);
    slope[0] = 0.03;
    slope[1] = 0.05;
    array_1d<double,3> height_gradient = -slope;

    WaveElementSteadyStateTest(manning, height, velocity, slope, height_gradient);
}

/**
 * @brief Check the WaveElement2D3N element still free surface and bottom topography x-gradient
 */
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3NBottomFriction, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.01;
    const double height = 5.0;
    array_1d<double,3> velocity = ZeroVector(3);
    velocity[0] = 0.6;
    array_1d<double,3> height_gradient = ZeroVector(3);
    height_gradient[0] = 0.03;
    const double central_height = height + height_gradient[0] / 3.; // at the barycenter of the element
    const array_1d<double,3> friction = std::pow(manning,2) * norm_2(velocity) * velocity / std::pow(central_height,4.0/3.0);
    const array_1d<double,3> slope = -height_gradient -friction;

    WaveElementSteadyStateTest(manning, height, velocity, slope, height_gradient);
}

} // namespace Testing

} // namespace Kratos
