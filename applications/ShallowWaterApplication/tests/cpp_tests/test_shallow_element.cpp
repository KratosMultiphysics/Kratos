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
#include "includes/properties.h"
#include "containers/model.h"
#include "custom_elements/shallow_element.h"
#include "shallow_water_application_variables.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType               IndexType;
typedef ModelPart::NodeIterator     NodeIteratorType;

/** 
 * Checks the ShallowWater2D3N element
 */
KRATOS_TEST_CASE_IN_SUITE(ShallowElement2D3N, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main");
    model_part.SetBufferSize(2);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(BATHYMETRY);
    model_part.AddNodalSolutionStepVariable(PROJECTED_VECTOR1);
    model_part.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
    model_part.AddNodalSolutionStepVariable(RAIN);

    // Process info creation
    const double delta_time = 0.1;
    const double dyn_tau = 0.005;
    const double gravity = 9.81;
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, dyn_tau);
    model_part.GetProcessInfo().SetValue(GRAVITY_Z, gravity);

    // Set the element properties
    Properties::Pointer property = model_part.pGetProperties(0);
    property->SetValue(MANNING, 0.004);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    model_part.CreateNewElement("ShallowElement2D3N", 1, elem_nodes, property);

    Element::Pointer element = model_part.pGetElement(1);

    array_1d<double,3> height;
    array_1d<double,3> bathymetry;
    height(0) = 7.0;
    height(1) = 9.0;
    height(2) = 6.0;
    bathymetry(0) = -7.0;
    bathymetry(1) = -8.0;
    bathymetry(2) = -6.0;
    // Set the nodal values
    for (IndexType i = 0; i < 3; i++)
    {
        element->GetGeometry()[i].FastGetSolutionStepValue(HEIGHT   ) = height(i);
        element->GetGeometry()[i].FastGetSolutionStepValue(HEIGHT, 1) = 0.0;
        element->GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_SCALAR1   ) = height(i) + 1;
        element->GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_SCALAR1, 1) = 0.0;
        element->GetGeometry()[i].FastGetSolutionStepValue(BATHYMETRY) = bathymetry(i);
        for (IndexType k = 0; k < 2; k++)
        {
            element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY   )[k] = 0.0;
            element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.0;
            element->GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_VECTOR1   )[k] = 0.0;
            element->GetGeometry()[i].FastGetSolutionStepValue(PROJECTED_VECTOR1, 1)[k] = 0.0;
        }
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(9);
    Matrix LHS = ZeroMatrix(9,9);

    element->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

    // Check the RHS values (the RHS is computed as the LHS x previous_solution, 
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    double tolerance = 1e-5;
    KRATOS_CHECK_NEAR(RHS(0),-1.63500, tolerance);
    KRATOS_CHECK_NEAR(RHS(1), 0.00000, tolerance);
    KRATOS_CHECK_NEAR(RHS(2), 1.66815, tolerance);
    KRATOS_CHECK_NEAR(RHS(3),-1.63500, tolerance);
    KRATOS_CHECK_NEAR(RHS(4), 0.00000, tolerance);
    KRATOS_CHECK_NEAR(RHS(5), 1.66518, tolerance);
    KRATOS_CHECK_NEAR(RHS(6),-1.63500, tolerance);
    KRATOS_CHECK_NEAR(RHS(7), 0.00000, tolerance);
    KRATOS_CHECK_NEAR(RHS(8), 1.66667, tolerance);
}

} // namespace Testing

} // namespace Kratos