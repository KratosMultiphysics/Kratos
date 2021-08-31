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
#include "includes/properties.h"
#include "custom_elements/swe.h"
#include "custom_elements/shallow_water_2d_3.h"
#include "shallow_water_application_variables.h"
#include "includes/mesh_moving_variables.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType               IndexType;
typedef ModelPart::NodeIterator     NodeIteratorType;

/**
 * Checks the ShallowWater2D3N element
 */
KRATOS_TEST_CASE_IN_SUITE(SWE2D3N, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MOMENTUM);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(FREE_SURFACE_ELEVATION);
    model_part.AddNodalSolutionStepVariable(TOPOGRAPHY);
    model_part.AddNodalSolutionStepVariable(PROJECTED_VECTOR1);
    model_part.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
    model_part.AddNodalSolutionStepVariable(RAIN);
    model_part.AddNodalSolutionStepVariable(MANNING);

    // Process info creation
    const double delta_time = 0.1;
    const double stab_factor = 0.005;
    const double gravity = 9.81;
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    model_part.GetProcessInfo().SetValue(STABILIZATION_FACTOR, stab_factor);
    model_part.GetProcessInfo().SetValue(GRAVITY_Z, gravity);
    model_part.GetProcessInfo().SetValue(DRY_HEIGHT, 0.1);
    model_part.GetProcessInfo().SetValue(PERMEABILITY, 0.1);
    model_part.GetProcessInfo().SetValue(DRY_DISCHARGE_PENALTY, 0.1);

    // Set the element properties
    Properties::Pointer property = model_part.CreateNewProperties(0);
    double manning = 0.004;
    property->SetValue(MANNING, manning);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    model_part.CreateNewElement("SWE2D3N", 1, elem_nodes, property);

    Element::Pointer element = model_part.pGetElement(1);

    array_1d<double,3> free_surface;
    array_1d<double,3> topography;
    free_surface(0) = 0.0;
    free_surface(1) = 1.0;
    free_surface(2) = 0.0;
    topography(0) = -7.0;
    topography(1) = -8.0;
    topography(2) = -6.0;
    // Set the nodal values
    for (IndexType i = 0; i < 3; i++)
    {
        element->GetGeometry()[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = free_surface(i);
        element->GetGeometry()[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION, 1) = free_surface(i) + 1;
        element->GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_SCALAR1   ) = free_surface(i) + 1;
        element->GetGeometry()[i].FastGetSolutionStepValue(TOPOGRAPHY) = topography(i);
        element->GetGeometry()[i].FastGetSolutionStepValue(MANNING) = manning;
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);

    const ProcessInfo& r_process_info = model_part.GetProcessInfo();
    element->CalculateLocalSystem(LHS, RHS, r_process_info);

    // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    double tolerance = 1e-5;
    KRATOS_CHECK_NEAR(RHS(0),-11.99000, tolerance);
    KRATOS_CHECK_NEAR(RHS(1),  0.00000, tolerance);
    KRATOS_CHECK_NEAR(RHS(2),  3.11838, tolerance);
    KRATOS_CHECK_NEAR(RHS(3),-11.99000, tolerance);
    KRATOS_CHECK_NEAR(RHS(4),  0.00000, tolerance);
    KRATOS_CHECK_NEAR(RHS(5),  0.21495, tolerance);
    KRATOS_CHECK_NEAR(RHS(6),-11.99000, tolerance);
    KRATOS_CHECK_NEAR(RHS(7),  0.00000, tolerance);
    KRATOS_CHECK_NEAR(RHS(8),  1.66667, tolerance);
}

void PerformSteadyStateTest(
    const double& rManning,
    const double& rHeight,
    const array_1d<double,3>& rMomentum,
    const array_1d<double,3>& rTopographySlope,
    const array_1d<double,3>& rHeightGradient = ZeroVector(3),
    const double& rTolerance = 1e-12)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MOMENTUM);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(TOPOGRAPHY);
    model_part.AddNodalSolutionStepVariable(RAIN);
    model_part.AddNodalSolutionStepVariable(MANNING);
    model_part.AddNodalSolutionStepVariable(WIND);
    model_part.AddNodalSolutionStepVariable(ATMOSPHERIC_PRESSURE);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(VERTICAL_VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_ACCELERATION);

    // Process info creation
    const double gravity = 9.81;
    const double stab_factor = 0.005;
    const double shock_stab_factor = 1.0;
    const double relative_dry_height = 0.1;
    const double density_water = 1000.0;
    const double density_air = 1.0;
    model_part.GetProcessInfo().SetValue(GRAVITY_Z, gravity);
    model_part.GetProcessInfo().SetValue(STABILIZATION_FACTOR, stab_factor);
    model_part.GetProcessInfo().SetValue(SHOCK_STABILIZATION_FACTOR, shock_stab_factor);
    model_part.GetProcessInfo().SetValue(RELATIVE_DRY_HEIGHT, relative_dry_height);
    model_part.GetProcessInfo().SetValue(DENSITY, density_water);
    model_part.GetProcessInfo().SetValue(DENSITY_AIR, density_air);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    Properties::Pointer property = model_part.CreateNewProperties(0);
    Element::Pointer element = model_part.CreateNewElement("ShallowWater2D3N", 1, elem_nodes, property);
    auto r_geom = element->GetGeometry();

    // Set the nodal values
    for (IndexType i = 0; i < 3; i++)
    {
        const array_1d<double,3> coords = r_geom[i].Coordinates();
        const auto height = rHeight + inner_prod(coords, rHeightGradient);
        const auto velocity = rMomentum / height;
        const auto topography = inner_prod(coords, rTopographySlope);

        r_geom[i].FastGetSolutionStepValue(HEIGHT) = height;
        r_geom[i].FastGetSolutionStepValue(VELOCITY) = velocity;
        r_geom[i].FastGetSolutionStepValue(MOMENTUM) = rMomentum;
        r_geom[i].FastGetSolutionStepValue(MANNING) = rManning;
        r_geom[i].FastGetSolutionStepValue(TOPOGRAPHY) = topography;
    }

    // Compute RHS and LHS
    Vector rhs = ZeroVector(9);
    Matrix lhs = ZeroMatrix(9,9);
    const ProcessInfo& r_process_info = model_part.GetProcessInfo();
    element->CalculateLocalSystem(lhs, rhs, r_process_info);

    // Check the RHS values. Since it is a steady solution the RHS must be zero
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(rhs, ZeroVector(9), rTolerance);
}

/**
 * @brief Check the ShallowWater2D3N element with still free surface
 */
KRATOS_TEST_CASE_IN_SUITE(SteadyStillSurfaceShallowWater2D3N, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    const array_1d<double,3> momentum = ZeroVector(3);
    const array_1d<double,3> slope = ZeroVector(3);

    PerformSteadyStateTest(manning, height, momentum, slope);
}

/**
 * @brief Check the ShallowWater2D3N element with subcritical x-aligned flow
 */
KRATOS_TEST_CASE_IN_SUITE(SteadySubcriticalShallowWater2D3N, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0328;
    const double height = 5.0;
    array_1d<double,3> momentum = ZeroVector(3);
    momentum[0] = 5.0;
    const array_1d<double,3> slope = -std::pow(manning, 2.) * norm_2(momentum) * momentum / std::pow(height, 10./3.);

    PerformSteadyStateTest(manning, height, momentum, slope);
}

/**
 * @brief Check the ShallowWater2D3N element with subcritical flow
 */
KRATOS_TEST_CASE_IN_SUITE(SteadySubcriticalSkewShallowWater2D3N, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0328;
    const double height = 5.0;
    array_1d<double,3> momentum = ZeroVector(3);
    momentum[0] = -3.0;
    momentum[1] = 4.0;
    const array_1d<double,3> slope = -std::pow(manning, 2.) * norm_2(momentum) * momentum / std::pow(height, 10./3.);

    PerformSteadyStateTest(manning, height, momentum, slope);
}

/**
 * @brief Check the ShallowWater2D3N element with subcritical flow (Froude = 2.72789)
 */
KRATOS_TEST_CASE_IN_SUITE(SteadySupercriticalShallowWater2D3N, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0328;
    const double height = 1.0;
    array_1d<double,3> momentum = ZeroVector(3);
    momentum[0] = 3.0;
    momentum[1] = 8.0;
    const array_1d<double,3> slope = -std::pow(manning, 2.) * norm_2(momentum) * momentum / std::pow(height, 10./3.);

    PerformSteadyStateTest(manning, height, momentum, slope);
}

/**
 * @brief Check the ShallowWater2D3N element with variable free surface
 * @detail A linear element can not pass this parcel test analytically, since the velocity is not linear.
 */
KRATOS_TEST_CASE_IN_SUITE(SteadyVariableFreeSurfaceShallowWater2D3N, ShallowWaterApplicationFastSuite)
{
    const double manning = 0.0;
    const double height = 5.0;
    array_1d<double,3> height_grad = ZeroVector(3);
    array_1d<double,3> momentum = ZeroVector(3);
    height_grad[0] = .05;
    momentum[0] = 5.0;
    momentum[1] = 0.0;
    const double gravity = 9.81;
    const double central_height = height + height_grad[0] / 3.; // at the barycenter of the element
    const array_1d<double,3> friction = std::pow(manning, 2.) * norm_2(momentum) * momentum / std::pow(central_height, 10./3.);
    const array_1d<double,3> slope = (inner_prod(momentum, momentum) / (gravity * std::pow(central_height, 3.)) -1.) * height_grad - friction;
    const double tolerance = 1e-3;

    PerformSteadyStateTest(manning, height, momentum, slope, height_grad, tolerance);
}

} // namespace Testing

} // namespace Kratos
