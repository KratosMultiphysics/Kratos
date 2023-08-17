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
    ShallowWaterTestsUtilities::AddVariables(model_part);

    // Process info creation
    ProcessInfo& r_process_info = model_part.GetProcessInfo();
    r_process_info.SetValue(GRAVITY_Z, 9.81);
    r_process_info.SetValue(STABILIZATION_FACTOR, 0.01);

    // Geometry creation
    ShallowWaterTestsUtilities::CreateGeometry(model_part, "WaveElement2D3N", "WaveCondition2D2N");

    // Set the nodal values
    for (auto& r_node : model_part.Nodes())
    {
        const array_1d<double,3> coords = r_node.Coordinates();
        const auto height = rHeight + inner_prod(coords, rHeightGradient);
        const auto topography = inner_prod(coords, rTopographySlope);

        r_node.FastGetSolutionStepValue(VELOCITY) = rVelocity;
        r_node.FastGetSolutionStepValue(HEIGHT) = height;
        r_node.FastGetSolutionStepValue(MANNING) = rManning;
        r_node.FastGetSolutionStepValue(TOPOGRAPHY) = topography;
    }

    // Compute RHS
    Vector rhs = ZeroVector(9);
    ShallowWaterTestsUtilities::CalculateAndAssembleRHS(model_part, rhs);

    // Check the RHS values. Since it is a steady solution the RHS must be zero
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(rhs, ZeroVector(9), rTolerance);
}

/**
 * @brief Check the WaveElement2D3N element with still free surface
 */
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3N_StillWater, ShallowWaterApplicationFastSuite)
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
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3N_StillTopographyGradient, ShallowWaterApplicationFastSuite)
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
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3N_StillTopographySkewGradient, ShallowWaterApplicationFastSuite)
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
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3N_VelocityAndGradient, ShallowWaterApplicationFastSuite)
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
KRATOS_TEST_CASE_IN_SUITE(WaveElement2D3N_BottomFriction, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.01;
    const double height = 5.0;
    array_1d<double,3> velocity = ZeroVector(3);
    velocity[0] = 0.2;
    array_1d<double,3> height_gradient = ZeroVector(3);
    height_gradient[0] = 0.01;
    const double central_height = height + height_gradient[0] / 3.; // at the barycenter of the element
    const array_1d<double,3> friction = std::pow(manning,2) * norm_2(velocity) * velocity / std::pow(central_height,4.0/3.0);
    const array_1d<double,3> slope = -height_gradient -friction;

    WaveElementSteadyStateTest(manning, height, velocity, slope, height_gradient);
}

} // namespace Testing

} // namespace Kratos
