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
#include "custom_elements/boussinesq_element.h"
#include "shallow_water_application_variables.h"

#include "utilities/math_utils.h"

namespace Kratos {

namespace Testing {

/**
 * Checks the BoussinesqElement2D3N
 */
KRATOS_TEST_CASE_IN_SUITE(BoussinesqElement2D3N, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(FREE_SURFACE_ELEVATION);
    model_part.AddNodalSolutionStepVariable(VERTICAL_VELOCITY);
    model_part.AddNodalSolutionStepVariable(VELOCITY_LAPLACIAN);
    model_part.AddNodalSolutionStepVariable(TOPOGRAPHY);
    model_part.AddNodalSolutionStepVariable(MANNING);

    // Process info creation
    const double delta_time = 0.1;
    const double stab_factor = 0.01;
    const double gravity = 9.81;
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    model_part.GetProcessInfo().SetValue(STABILIZATION_FACTOR, stab_factor);
    model_part.GetProcessInfo().SetValue(GRAVITY_Z, gravity);
    model_part.GetProcessInfo().SetValue(DRY_HEIGHT, 0.1);

    // Set the element properties
    Properties::Pointer property = model_part.CreateNewProperties(0);
    double manning = 0.004;
    property->SetValue(MANNING, manning);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<std::size_t> elem_nodes {1, 2, 3};
    model_part.CreateNewElement("BoussinesqElement2D3N", 1, elem_nodes, property);

    Element::Pointer element = model_part.pGetElement(1);

    array_1d<double,3> topography;
    array_1d<double,3> free_surface;
    array_1d<array_1d<double,3>,3> velocity;
    array_1d<array_1d<double,3>,3> velocity_laplacian;
    array_1d<array_1d<double,3>,3> acceleration;
    array_1d<double,3> vertical_velocity;
    topography[0] = -1.0;
    topography[1] = -1.0;
    topography[2] = -1.0;
    free_surface[0] = 0.08;
    free_surface[1] = 0.06;
    free_surface[2] = 0.08;
    velocity[0] = array_1d<double,3>({0.2, 0.0, 0.0});
    velocity[1] = array_1d<double,3>({0.1, 0.0, 0.0});
    velocity[2] = array_1d<double,3>({0.2, 0.0, 0.0});
    velocity_laplacian[0] = array_1d<double,3>({0.0, 0.0, 0.0});
    velocity_laplacian[1] = array_1d<double,3>({0.05, 0.0, 0.0});
    velocity_laplacian[2] = array_1d<double,3>({0.0, 0.0, 0.0});
    acceleration[0] = array_1d<double,3>({0.212862, 0.0, 0.0});
    acceleration[1] = array_1d<double,3>({0.212867, 0.0, 0.0});
    acceleration[2] = array_1d<double,3>({0.212861, 0.0, 0.0});
    vertical_velocity[0] = 0.113501;
    vertical_velocity[1] = 0.113501;
    vertical_velocity[2] = 0.113501;
    // Set the nodal values
    for (std::size_t i = 0; i < 3; i++)
    {
        element->GetGeometry()[i].FastGetSolutionStepValue(TOPOGRAPHY) = topography[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = free_surface[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY) = velocity[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN) = velocity_laplacian[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION) = acceleration[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(VERTICAL_VELOCITY) = vertical_velocity[i];
    }

    // Calculate the local system, M*dot(u) = F(u)
    Matrix M;
    Matrix LHS;
    Vector RHS;
    Vector derivatives;
    const ProcessInfo& r_process_info = model_part.GetProcessInfo();

    element->CalculateMassMatrix(M, r_process_info);
    element->CalculateLocalSystem(LHS, RHS, r_process_info);
    element->GetFirstDerivativesVector(derivatives);

    double tolerance = 1e-6;
    Vector increment = prod(M,derivatives) -RHS;
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(increment, ZeroVector(9), tolerance);
}

} // namespace Testing

} // namespace Kratos
