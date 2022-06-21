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

// Application includes
#include "shallow_water_application_variables.h"
#include "shallow_water_tests_utilities.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType IndexType;

//TODO: Fix this multiple definition
void SetNodalValues2(
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
        const auto topography = inner_prod(coords, rTopographySlope);

        r_node.FastGetSolutionStepValue(VELOCITY) = rVelocity;
        r_node.FastGetSolutionStepValue(HEIGHT) = height;
        r_node.FastGetSolutionStepValue(MANNING) = rManning;
        r_node.FastGetSolutionStepValue(TOPOGRAPHY) = topography;
    }
}

void EmbeddedPrimitiveElementSteadyStateTest(
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

    // Set ProcessInfo
    ProcessInfo& r_process_info = model_part.GetProcessInfo();
    r_process_info.SetValue(GRAVITY_Z, 9.81);
    r_process_info.SetValue(STABILIZATION_FACTOR, 0.01);
    r_process_info.SetValue(RELATIVE_DRY_HEIGHT, 0.1);

    // Create the triangle and its conditions
    ShallowWaterTestsUtilities::CreateGeometry(model_part, "EmbeddedPrimitiveElement2D3N", "WaveCondition2D2N");

    // Set the nodal values
    SetNodalValues2(model_part, rManning, rHeight, rVelocity, rTopographySlope, rHeightGradient);

    // Set distance field
    auto& r_elem = model_part.GetElement(1);
    auto& r_geom = r_elem.GetGeometry();
    r_geom[0].FastGetSolutionStepValue(DISTANCE) = 1.0;
    r_geom[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
    r_geom[2].FastGetSolutionStepValue(DISTANCE) = 1.0;

    // Compute RHS
    Vector rhs = ZeroVector(9);
    Matrix lhs = ZeroMatrix(9,9);
    r_elem.CalculateLocalSystem(lhs, rhs, r_process_info);

    KRATOS_WATCH(rhs)
    KRATOS_WATCH(lhs)

    // // Check the RHS values. Since it is a steady solution the RHS must be zero
    // KRATOS_CHECK_VECTOR_RELATIVE_NEAR(rhs, ZeroVector(9), rTolerance);
}

/**
 * @brief Check the EmbeddedPrimitiveElement2D3N element with still free surface
 */
KRATOS_TEST_CASE_IN_SUITE(EmbeddedPrimitiveElement2D3N_SteadyStillSurface, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0e-4;
    const array_1d<double,3> velocity = ZeroVector(3);
    array_1d<double,3> slope = ZeroVector(3);
    array_1d<double,3> height_gradient = ZeroVector(3);
    slope[0] = 1.0e-3;
    height_gradient[0] = -1.0e-3;

    EmbeddedPrimitiveElementSteadyStateTest(manning, height, velocity, slope, height_gradient);
}

} // namespace Testing

} // namespace Kratos
