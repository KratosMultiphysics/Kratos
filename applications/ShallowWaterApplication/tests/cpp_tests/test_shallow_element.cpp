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
#include "includes/model_part.h"
#include "custom_elements/shallow_element.h"
#include "shallow_water_application_variables.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::IndexType               IndexType;
typedef ModelPart::NodeIterator     NodeIteratorType;

/** 
 * Checks the ShallowWater2D3N element
 */
KRATOS_TEST_CASE_IN_SUITE(ShallowElement2D34, ShallowWaterApplicationFastSuite)
{
    ModelPart model_part("main");
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

    Element::Pointer element = modelPart.pGetElement(1);

}

} // namespace Testing

} // namespace Kratos