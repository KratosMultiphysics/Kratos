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
#include "custom_utilities/shallow_water_tests_utilities.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType IndexType;

void ConservativeElementSteadyStateTest(
    const double& rManning,
    const double& rHeight,
    const array_1d<double,3>& rMomentum,
    const array_1d<double,3>& rTopographySlope,
    const array_1d<double,3>& rHeightGradient = ZeroVector(3),
    const double& rTolerance = 1e-12)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);

    // Variables
    ShallowWaterTestsUtilities::AddVariables(model_part);

    // Set ProcessInfo
    ProcessInfo& r_process_info = model_part.GetProcessInfo();
    r_process_info.SetValue(GRAVITY_Z, 9.81);
    r_process_info.SetValue(STABILIZATION_FACTOR, 0.01);
    r_process_info.SetValue(RELATIVE_DRY_HEIGHT, 0.1);

    // Geometry creation
    // TestCreateGeometry(model_part, "ConservativeElement2D3N", "ConservativeCondition2D3N");
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Properties::Pointer property = model_part.CreateNewProperties(0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    Element::Pointer element = model_part.CreateNewElement("ConservativeElement2D3N", 1, elem_nodes, property);
    std::vector<IndexType> cond_1_nodes {1, 2};
    std::vector<IndexType> cond_2_nodes {2, 3};
    std::vector<IndexType> cond_3_nodes {3, 1};
    Condition::Pointer cond_1 = model_part.CreateNewCondition("ConservativeCondition2D2N", 1, cond_1_nodes, property);
    Condition::Pointer cond_2 = model_part.CreateNewCondition("ConservativeCondition2D2N", 2, cond_2_nodes, property);
    Condition::Pointer cond_3 = model_part.CreateNewCondition("ConservativeCondition2D2N", 3, cond_3_nodes, property);

    // Set the nodal values
    auto r_geom = element->GetGeometry();
    for (auto& r_node : r_geom)
    {
        const array_1d<double,3> coords = r_node.Coordinates();
        const auto height = rHeight + inner_prod(coords, rHeightGradient);
        const auto velocity = rMomentum / height;
        const auto topography = inner_prod(coords, rTopographySlope);

        r_node.FastGetSolutionStepValue(MOMENTUM) = rMomentum;
        r_node.FastGetSolutionStepValue(VELOCITY) = velocity;
        r_node.FastGetSolutionStepValue(HEIGHT) = height;
        r_node.FastGetSolutionStepValue(MANNING) = rManning;
        r_node.FastGetSolutionStepValue(TOPOGRAPHY) = topography;
    }

    // Compute RHS
    Vector rhs = ZeroVector(9);
    Matrix lhs = ZeroMatrix(9,9);
    element->CalculateLocalSystem(lhs, rhs, r_process_info);

    Vector rhs_cond = ZeroVector(6);
    cond_1->CalculateRightHandSide(rhs_cond, r_process_info);
    ShallowWaterTestsUtilities::AssembleRHS(rhs, rhs_cond, cond_1_nodes);

    cond_2->CalculateRightHandSide(rhs_cond, r_process_info);
    ShallowWaterTestsUtilities::AssembleRHS(rhs, rhs_cond, cond_2_nodes);

    cond_3->CalculateRightHandSide(rhs_cond, r_process_info);
    ShallowWaterTestsUtilities::AssembleRHS(rhs, rhs_cond, cond_3_nodes);

    // Check the RHS values. Since it is a steady solution the RHS must be zero
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(rhs, ZeroVector(9), rTolerance);
}

/**
 * @brief Check the ConservativeElement2D3N element with still free surface
 */
KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3N_SteadyStillSurface, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    const array_1d<double,3> momentum = ZeroVector(3);
    const array_1d<double,3> slope = ZeroVector(3);

    ConservativeElementSteadyStateTest(manning, height, momentum, slope);
}

/**
 * @brief Check the ConservativeElement2D3N element still free surface and bottom topography x-gradient
 */
KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3N_SteadyTopographyGradient, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    array_1d<double,3> momentum = ZeroVector(3);
    array_1d<double,3> slope = ZeroVector(3);
    slope[0] = 0.05;
    array_1d<double,3> height_gradient = -slope;

    ConservativeElementSteadyStateTest(manning, height, momentum, slope, height_gradient);
}

/**
 * @brief Check the ConservativeElement2D3N element still free surface and bottom topography skew gradient
 */
KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3N_SteadyTopographySkewGradient, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    array_1d<double,3> momentum = ZeroVector(3);
    array_1d<double,3> slope = ZeroVector(3);
    slope[0] = 0.05;
    slope[1] = 0.05;
    array_1d<double,3> height_gradient = -slope;

    ConservativeElementSteadyStateTest(manning, height, momentum, slope, height_gradient);
}

/**
 * @brief Check the ConservativeElement2D3N element steady subcritical flow
 */
KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3N_SteadySubcriticalFlow, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0328;
    const double height = 5.0;
    array_1d<double,3> momentum = ZeroVector(3);
    momentum[0] = -3.0;//0.6;
    momentum[1] = 4.0;//1.0;
    const array_1d<double,3> slope = -std::pow(manning, 2.) * norm_2(momentum) * momentum / std::pow(height, 10./3.);

    ConservativeElementSteadyStateTest(manning, height, momentum, slope);
}

/**
 * @brief Check the ConservativeElement2D3N element steady state with variable velocity
 */
KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3N_SteadyStateVariableVelocityX, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.01;
    const double height = 5.0;
    array_1d<double,3> momentum = ZeroVector(3);
    momentum[0] = 0.6;
    array_1d<double,3> height_grad = ZeroVector(3);
    height_grad[0] = 0.03;
    const double gravity = 9.81;
    const double central_height = height + height_grad[0] / 3.; // at the barycenter of the element
    const array_1d<double,3> friction = std::pow(manning,2) * norm_2(momentum) * momentum / std::pow(central_height,10.0/3.0);
    const array_1d<double,3> slope = (inner_prod(momentum, momentum) / (gravity * std::pow(central_height, 3.)) -1.) * height_grad -friction;
    const double tolerance = 1e-8;

    ConservativeElementSteadyStateTest(manning, height, momentum, slope, height_grad, tolerance);
}

/**
 * @brief Check the ConservativeElement2D3N element steady state with variable velocity
 */
KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3N_SteadyStateVariableVelocityY, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.01;
    const double height = 5.0;
    array_1d<double,3> momentum = ZeroVector(3);
    momentum[1] = 0.2;
    array_1d<double,3> height_grad = ZeroVector(3);
    height_grad[1] = 0.01;
    const double gravity = 9.81;
    const double central_height = height + height_grad[0] / 3.; // at the barycenter of the element
    const array_1d<double,3> friction = std::pow(manning,2) * norm_2(momentum) * momentum / std::pow(central_height,10.0/3.0);
    const array_1d<double,3> slope = (inner_prod(momentum, momentum) / (gravity * std::pow(central_height, 3.)) -1.) * height_grad -friction;
    const double tolerance = 1e-8;

    ConservativeElementSteadyStateTest(manning, height, momentum, slope, height_grad, tolerance);
}

} // namespace Testing

} // namespace Kratos
