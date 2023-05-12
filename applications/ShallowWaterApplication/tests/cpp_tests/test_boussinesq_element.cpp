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
// #include "utilities/math_utils.h"

namespace Kratos {

namespace Testing {

void InitializeModelPart(ModelPart& rModelPart)
{
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(MOMENTUM);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(HEIGHT);
    rModelPart.AddNodalSolutionStepVariable(FREE_SURFACE_ELEVATION);
    rModelPart.AddNodalSolutionStepVariable(VERTICAL_VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(DISPERSION_H);
    rModelPart.AddNodalSolutionStepVariable(DISPERSION_V);
    rModelPart.AddNodalSolutionStepVariable(TOPOGRAPHY);
    rModelPart.AddNodalSolutionStepVariable(MANNING);

    // Process info creation
    const double delta_time = 0.1;
    const double stab_factor = 0.01;
    const double gravity = 9.81;
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, delta_time);
    rModelPart.GetProcessInfo().SetValue(STABILIZATION_FACTOR, stab_factor);
    rModelPart.GetProcessInfo().SetValue(GRAVITY_Z, gravity);
    rModelPart.GetProcessInfo().SetValue(DRY_HEIGHT, 0.1);

    // Set the element properties
    Properties::Pointer property = rModelPart.CreateNewProperties(0);
    double manning = 0.004;
    property->SetValue(MANNING, manning);
}

/**
 * Check the BoussinesqElement2D3N
 */
KRATOS_TEST_CASE_IN_SUITE(BoussinesqElement2D3N_FlatBottom, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);

    InitializeModelPart(model_part);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<std::size_t> elem_nodes {1, 2, 3};
    model_part.CreateNewElement("BoussinesqElement2D3N", 1, elem_nodes, model_part.pGetProperties(0));

    Element::Pointer element = model_part.pGetElement(1);

    constexpr std::size_t TNumNodes = 3;
    array_1d<double,TNumNodes> topography;
    array_1d<double,TNumNodes> free_surface;
    array_1d<array_1d<double,3>,TNumNodes> velocity;
    array_1d<array_1d<double,3>,TNumNodes> dispersion;
    array_1d<array_1d<double,3>,TNumNodes> acceleration;
    array_1d<double,TNumNodes> vertical_velocity;
    topography[0] = -1.0;
    topography[1] = -1.0;
    topography[2] = -1.0;
    free_surface[0] = 0.08;
    free_surface[1] = 0.06;
    free_surface[2] = 0.08;
    velocity[0] = array_1d<double,3>({0.2, 0.0, 0.0});
    velocity[1] = array_1d<double,3>({0.1, 0.0, 0.0});
    velocity[2] = array_1d<double,3>({0.2, 0.0, 0.0});
    dispersion[0] = array_1d<double,3>({0.0, 0.0, 0.0});
    dispersion[1] = array_1d<double,3>({0.01, 0.0, 0.0});
    dispersion[2] = array_1d<double,3>({0.0, 0.0, 0.0});
    acceleration[0] = array_1d<double,3>({0.213822,-0.00237072, 0.0});
    acceleration[1] = array_1d<double,3>({0.208574, 0.00000000, 0.0});
    acceleration[2] = array_1d<double,3>({0.216193, 0.00237072, 0.0});
    vertical_velocity[0] = 0.102002;
    vertical_velocity[1] = 0.097998;
    vertical_velocity[2] = 0.102002;

    // Set the nodal values
    for (std::size_t i = 0; i < element->GetGeometry().size(); i++)
    {
        element->GetGeometry()[i].FastGetSolutionStepValue(TOPOGRAPHY) = topography[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(HEIGHT) = free_surface[i] - topography[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY) = velocity[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(DISPERSION_H) = dispersion[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(DISPERSION_V) = dispersion[i];
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

    // Matrix inv_M;
    // double det_M;
    // MathUtils<double>::InvertMatrix(M, inv_M, det_M);
    // KRATOS_WATCH_CERR(prod(inv_M, RHS))

    double tolerance = 1e-6;
    Vector increment = prod(M,derivatives) -RHS;
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(increment, ZeroVector(9), tolerance);
}

/**
 * Check the BoussinesqElement2D4N
 */
KRATOS_TEST_CASE_IN_SUITE(BoussinesqElement2D4N_FlatBottom, ShallowWaterApplicationFastSuite)
{
    Model model;
    ModelPart& model_part = model.CreateModelPart("main", 2);

    InitializeModelPart(model_part);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
    std::vector<std::size_t> elem_nodes {1, 2, 3, 4};
    model_part.CreateNewElement("BoussinesqElement2D4N", 1, elem_nodes, model_part.pGetProperties(0));

    Element::Pointer element = model_part.pGetElement(1);

    constexpr std::size_t TNumNodes = 4;
    array_1d<double,TNumNodes> topography;
    array_1d<double,TNumNodes> free_surface;
    array_1d<array_1d<double,3>,TNumNodes> velocity;
    array_1d<array_1d<double,3>,TNumNodes> dispersion;
    array_1d<array_1d<double,3>,TNumNodes> acceleration;
    array_1d<double,TNumNodes> vertical_velocity;
    topography[0] = -1.0;
    topography[1] = -1.0;
    topography[2] = -1.0;
    topography[3] = -1.0;
    free_surface[0] = 0.08;
    free_surface[1] = 0.08;
    free_surface[2] = 0.06;
    free_surface[3] = 0.06;
    velocity[0] = array_1d<double,3>({0.0, 0.2, 0.0});
    velocity[1] = array_1d<double,3>({0.0, 0.2, 0.0});
    velocity[2] = array_1d<double,3>({0.0, 0.1, 0.0});
    velocity[3] = array_1d<double,3>({0.0, 0.1, 0.0});
    dispersion[0] = array_1d<double,3>({0.0, 0.0, 0.0});
    dispersion[1] = array_1d<double,3>({0.0, 0.0, 0.0});
    dispersion[2] = array_1d<double,3>({0.0, 0.05, 0.0});
    dispersion[3] = array_1d<double,3>({0.0, 0.05, 0.0});
    acceleration[0] = array_1d<double,3>({-0.002258, 0.213936, 0.0});
    acceleration[1] = array_1d<double,3>({ 0.002258, 0.213936, 0.0});
    acceleration[2] = array_1d<double,3>({ 0.002256, 0.208458, 0.0});
    acceleration[3] = array_1d<double,3>({-0.002256, 0.208458, 0.0});
    vertical_velocity[0] = 0.0620025;
    vertical_velocity[1] = 0.0620025;
    vertical_velocity[2] = 0.0579975;
    vertical_velocity[3] = 0.0579975;

    // Set the nodal values
    for (std::size_t i = 0; i < element->GetGeometry().size(); i++)
    {
        element->GetGeometry()[i].FastGetSolutionStepValue(TOPOGRAPHY) = topography[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(HEIGHT) = free_surface[i] - topography[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY) = velocity[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(DISPERSION_H) = dispersion[i];
        element->GetGeometry()[i].FastGetSolutionStepValue(DISPERSION_V) = dispersion[i];
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

    // Matrix inv_M;
    // double det_M;
    // MathUtils<double>::InvertMatrix(M, inv_M, det_M);
    // KRATOS_WATCH_CERR(prod(inv_M, RHS))

    double tolerance = 1e-6;
    Vector increment = prod(M,derivatives) -RHS;
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(increment, ZeroVector(12), tolerance);
}

} // namespace Testing

} // namespace Kratos
