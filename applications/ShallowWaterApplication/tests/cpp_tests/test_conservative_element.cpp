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
#include "shallow_water_tests_utilities.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType IndexType;

void SetNodalValues(
    ModelPart& rModelPart,
    const double& rManning,
    const double& rHeight,
    const array_1d<double,3>& rMomentum,
    const array_1d<double,3>& rTopographySlope,
    const array_1d<double,3>& rHeightGradient)
{
    for (auto& r_node : rModelPart.Nodes())
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
}

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

    // Create the triangle and its conditions
    ShallowWaterTestsUtilities::CreateGeometry(model_part, "ConservativeElementRV2D3N", "ConservativeCondition2D2N");

    // Set the nodal values
    SetNodalValues(model_part, rManning, rHeight, rMomentum, rTopographySlope, rHeightGradient);

    // Compute RHS
    Vector rhs = ZeroVector(9);
    ShallowWaterTestsUtilities::CalculateAndAssembleRHS(model_part, rhs);

    // Check the RHS values. Since it is a steady solution the RHS must be zero
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(rhs, ZeroVector(9), rTolerance);
}

void ConservativeElementSteadyStateTestParts(
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
    r_process_info.SetValue(INTEGRATE_BY_PARTS, true);
    r_process_info.SetValue(GRAVITY_Z, 9.81);
    r_process_info.SetValue(STABILIZATION_FACTOR, 0.01);
    r_process_info.SetValue(RELATIVE_DRY_HEIGHT, 0.1);

    // Create the triangle and its conditions
    ShallowWaterTestsUtilities::CreateGeometry(model_part, "ConservativeElementRV2D3N", "ConservativeCondition2D2N");

    // Set the nodal values
    SetNodalValues(model_part, rManning, rHeight, rMomentum, rTopographySlope, rHeightGradient);

    // Compute RHS
    Vector rhs = ZeroVector(9);
    ShallowWaterTestsUtilities::CalculateAndAssembleRHS(model_part, rhs);

    // Check the RHS values. Since it is a steady solution the RHS must be zero
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(rhs, ZeroVector(9), rTolerance);
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
    momentum[0] = 0.2;
    array_1d<double,3> height_grad = ZeroVector(3);
    height_grad[0] = 0.01;
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

/**
 * @brief Check the ConservativeElement2D3N element integrated by parts with still free surface
 */
KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3NByParts_SteadyStillSurface, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    const array_1d<double,3> momentum = ZeroVector(3);
    const array_1d<double,3> slope = ZeroVector(3);

    ConservativeElementSteadyStateTestParts(manning, height, momentum, slope);
}

/**
 * @brief Check the ConservativeElement2D3N element integrated by parts with still free surface and bottom topography skew gradient
 */
// KRATOS_TEST_CASE_IN_SUITE(ConservativeElement2D3NByParts_SteadyTopographySkewGradient, ShallowWaterApplicationFastSuite)
// {
//     const double manning = 0.0;
//     const double height = 5.0;
//     array_1d<double,3> momentum = ZeroVector(3);
//     array_1d<double,3> slope = ZeroVector(3);
//     slope[0] = 0.03;
//     slope[1] = 0.03;
//     array_1d<double,3> height_gradient = -slope;

//     ConservativeElementSteadyStateTestParts(manning, height, momentum, slope, height_gradient);
// }

} // namespace Testing

} // namespace Kratos
