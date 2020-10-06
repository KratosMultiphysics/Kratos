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
#include "shallow_water_application_variables.h"

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
    model_part.AddNodalSolutionStepVariable(EQUIVALENT_MANNING);
    model_part.AddNodalSolutionStepVariable(POROSITY);

    // Process info creation
    const double delta_time = 0.1;
    const double stab_factor = 0.005;
    const double gravity = 9.81;
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    model_part.GetProcessInfo().SetValue(STABILIZATION_FACTOR, stab_factor);
    model_part.GetProcessInfo().SetValue(GRAVITY_Z, gravity);
    model_part.GetProcessInfo().SetValue(DRY_HEIGHT, 0.1);
    model_part.GetProcessInfo().SetValue(WATER_HEIGHT_UNIT_CONVERTER, 1.0);
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
        element->GetGeometry()[i].FastGetSolutionStepValue(POROSITY) = 1.0;
        element->GetGeometry()[i].FastGetSolutionStepValue(EQUIVALENT_MANNING) = manning;
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);

    element->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

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

} // namespace Testing

} // namespace Kratos
