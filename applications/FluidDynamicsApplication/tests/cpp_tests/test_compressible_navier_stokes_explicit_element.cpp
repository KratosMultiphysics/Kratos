//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "containers/model.h"
#include "custom_elements/compressible_navier_stokes_explicit.h"

namespace Kratos {
namespace Testing {

/**
 * @brief Create a Compressible Navier Stokes Explicit 2D3N Element object
 * This method creates a compressible Navier-Stokes triangular element to be tested
 * @param rModelPart Referente to the model part in which the element is created
 */
void CreateCompressibleNavierStokesExplicit2D3NElement(ModelPart& rModelPart)
{
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(REACTION_DENSITY);
    rModelPart.AddNodalSolutionStepVariable(MOMENTUM);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(REACTION_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    rModelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);
    rModelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);

    // Set the element properties
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
    p_elem_prop->SetValue(CONDUCTIVITY,  0.0257);
    p_elem_prop->SetValue(SPECIFIC_HEAT,  718);
    p_elem_prop->SetValue(HEAT_CAPACITY_RATIO, 1.4);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    auto p_elem_expl = rModelPart.CreateNewElement("CompressibleNavierStokesExplicit2D3N", 1, elem_nodes, p_elem_prop);
}

/**
 * @brief Create a Compressible Navier Stokes Explicit 3D4N Element object
 * This method creates a compressible Navier-Stokes tetrahedron element to be tested
 * @param rModelPart Referente to the model part in which the element is created
 */
void CreateCompressibleNavierStokesExplicit3D4NElement(ModelPart& rModelPart)
{
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(BODY_FORCE);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(REACTION_DENSITY);
    rModelPart.AddNodalSolutionStepVariable(MOMENTUM);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(REACTION_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    rModelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);
    rModelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);

    // Set the element properties
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
    p_elem_prop->SetValue(CONDUCTIVITY,  0.0257);
    p_elem_prop->SetValue(SPECIFIC_HEAT,  718);
    p_elem_prop->SetValue(HEAT_CAPACITY_RATIO, 1.4);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3, 4};
    auto p_elem_expl = rModelPart.CreateNewElement("CompressibleNavierStokesExplicit3D4N", 1, elem_nodes, p_elem_prop);
}

/**
 * @brief Test the 2D explicit compressible Navier-Stokes element RHS
 * This test compares the 2D explicit RHS contribution to the implicit one without the inertial terms
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicitVsImplicitRHS2D3N, FluidDynamicsApplicationFastSuite)
{
    // Create the test geometry
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 2);
    CreateCompressibleNavierStokesExplicit2D3NElement(r_model_part);

    // Define and set the nodal values
    double R = 287.05;						// J/(kg*K)
    double T = 0.00247;                     // Temperature in K
    array_1d<double, 3> momentum;
    array_1d<double, 3> f_ext;
    f_ext[0] = 0.0; // x-volume force
    f_ext[1] = 0.0; // y-volume force
    f_ext[2] = 0.0; // z-volume force
    double r = 0.0; // external pressure

    for (auto &r_node : r_model_part.Nodes()){
        const double velocity = 2.9 * (1.0 - (r_node.Id() / 10.0));
        const double density = 1.772 * (1.0 - (r_node.Id() / 10.0));
        momentum[0] = velocity * density;
        momentum[1] = velocity * density;
        momentum[2] = 0.0;
        const double total_energy = (density*R*T/(1.4-1)) + ((inner_prod(momentum,momentum))/(2*density));
        r_node.FastGetSolutionStepValue(DENSITY) = density;
        r_node.SetValue(DENSITY_TIME_DERIVATIVE, 0.0);
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.SetValue(MOMENTUM_TIME_DERIVATIVE, MOMENTUM_TIME_DERIVATIVE.Zero());
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        r_node.SetValue(TOTAL_ENERGY_TIME_DERIVATIVE, 0.0);
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
    auto p_elem_expl = r_model_part.pGetElement(1);
    p_elem_expl->Initialize(r_process_info);
    p_elem_expl->AddExplicitContribution(r_process_info);

    // Check obtained RHS values
    // We check against the RHS obtained with the original implicit element without the intertial terms and the shock capturing
    std::vector<double> RHS_impl({0.134831, 0.328435, 0.324492, 1.15847, -0.695681, -2.39476, -2.34955, -8.87494, -0.698156, -2.30026, -2.39888, -8.75893});
    std::vector<double> RHS_expl(12);
    for (unsigned int i_node = 0; i_node < r_model_part.NumberOfNodes(); ++i_node) {
        const auto it_node = r_model_part.NodesBegin() + i_node;
        RHS_expl[i_node * 4    ] = it_node->FastGetSolutionStepValue(REACTION_DENSITY);
        RHS_expl[i_node * 4 + 1] = it_node->FastGetSolutionStepValue(REACTION_X);
        RHS_expl[i_node * 4 + 2] = it_node->FastGetSolutionStepValue(REACTION_Y);
        RHS_expl[i_node * 4 + 3] = it_node->FastGetSolutionStepValue(REACTION_ENERGY);
    }
    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_impl, 1e-5);
}

/**
 * @brief Test the 2D explicit compressible Navier-Stokes element RHS
 * This test checks the residual calculation of the explicit compressible Navier-Stokes triangular element
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicitRHS2D3N, FluidDynamicsApplicationFastSuite)
{
    // Create the test geometry
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 2);
    CreateCompressibleNavierStokesExplicit2D3NElement(r_model_part);

    // Define and set the nodal values
    const double R = 287.05; // J/(kg*K)
    const double T = 0.00247; // Temperature in K
    array_1d<double, 3> momentum;
    array_1d<double, 3> f_ext;
    f_ext[0] = 0.0; // x-volume force
    f_ext[1] = 0.0; // y-volume force
    f_ext[2] = 0.0; // z-volume force
    const double r = 0.0; // external pressure
    const double dt = 1.0e-1; // fake delta time to calculate the time derivatives

    for (auto &r_node : r_model_part.Nodes()){
        const double velocity = 2.9 * (1.0 - (r_node.Id() / 10.0));
        const double density = 1.772 * (1.0 - (r_node.Id() / 10.0));
        momentum[0] = velocity * density;
        momentum[1] = velocity * density;
        momentum[2] = 0.0;
        const double total_energy = (density*R*T/(1.4-1)) + ((inner_prod(momentum,momentum))/(2*density));
        r_node.FastGetSolutionStepValue(DENSITY) = density;
        r_node.SetValue(DENSITY_TIME_DERIVATIVE, density / dt);
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.SetValue(MOMENTUM_TIME_DERIVATIVE, momentum / dt);
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        r_node.SetValue(TOTAL_ENERGY_TIME_DERIVATIVE, total_energy / dt);
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
    auto p_elem_expl = r_model_part.pGetElement(1);
    p_elem_expl->Initialize(r_process_info);
    p_elem_expl->AddExplicitContribution(r_process_info);

    // Check obtained RHS values
    // We check against the RHS obtained with the original implicit element without the intertial terms and the shock capturing
    std::vector<double> RHS_ref({-1.88987, -4.7346, -4.73855, -15.1509, 0.316671, 0.29379, 0.0249656, -0.720196, 0.314195, 0.0742576, 0.289669, -0.604187});
    std::vector<double> RHS_expl(12);
    for (unsigned int i_node = 0; i_node < r_model_part.NumberOfNodes(); ++i_node) {
        const auto it_node = r_model_part.NodesBegin() + i_node;
        RHS_expl[i_node * 4    ] = it_node->FastGetSolutionStepValue(REACTION_DENSITY);
        RHS_expl[i_node * 4 + 1] = it_node->FastGetSolutionStepValue(REACTION_X);
        RHS_expl[i_node * 4 + 2] = it_node->FastGetSolutionStepValue(REACTION_Y);
        RHS_expl[i_node * 4 + 3] = it_node->FastGetSolutionStepValue(REACTION_ENERGY);
    }
    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_ref, 1e-4);
}

/**
 * @brief Test the 3D explicit compressible Navier-Stokes element RHS
 * This test compares the 3D explicit RHS contribution to the implicit one without the inertial terms
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicitVsImplicitRHS3D4N, FluidDynamicsApplicationFastSuite)
{
    // Supersonic test for constant variables
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 2);
    CreateCompressibleNavierStokesExplicit3D4NElement(r_model_part);

    // Define and set the nodal values
    double R = 287.05;						// J/(kg*K)
    double T = 0.00247;                     // Temperature in K
    array_1d<double, 3> momentum;
    array_1d<double, 3> f_ext;
    f_ext[0] = 0.0; // x-volume force
    f_ext[1] = 0.0; // y-volume force
    f_ext[2] = 0.0; // z-volume force
    double r = 0.0; // external pressure

    for (auto &r_node : r_model_part.Nodes()){
        const double velocity = 2.9 * (1.0 - (r_node.Id() / 10.0));
        const double density = 1.772 * (1.0 - (r_node.Id() / 10.0));
        momentum[0] = velocity * density;
        momentum[1] = velocity * density;
        momentum[2] = velocity * density;
        const double total_energy = (density*R*T/(1.4-1)) + ((inner_prod(momentum,momentum))/(2*density));
        r_node.FastGetSolutionStepValue(DENSITY) = density;
        r_node.SetValue(DENSITY_TIME_DERIVATIVE, 0.0);
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.SetValue(MOMENTUM_TIME_DERIVATIVE, MOMENTUM_TIME_DERIVATIVE.Zero());
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        r_node.SetValue(TOTAL_ENERGY_TIME_DERIVATIVE, 0.0);
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
    auto p_elem_expl= r_model_part.pGetElement(1);
    p_elem_expl->Initialize(r_process_info);
    p_elem_expl->AddExplicitContribution(r_process_info);

    // Check obtained RHS values
    // We check against the RHS obtained with the original implicit element without the intertial terms and the shock capturing
    std::vector<double> RHS_impl({0.186062, 0.525773, 0.525316, 0.524028, 2.61985, -0.329746, -1.07454, -1.04549, -1.04842, -5.41058, -0.33065, -1.02446, -1.07384, -1.03918, -5.33134, -0.330744, -1.01243, -1.02437, -1.06775, -5.25722});
    std::vector<double> RHS_expl(20);
    for (unsigned int i_node = 0; i_node < r_model_part.NumberOfNodes(); ++i_node) {
        const auto it_node = r_model_part.NodesBegin() + i_node;
        RHS_expl[i_node * 5    ] = it_node->FastGetSolutionStepValue(REACTION_DENSITY);
        RHS_expl[i_node * 5 + 1] = it_node->FastGetSolutionStepValue(REACTION_X);
        RHS_expl[i_node * 5 + 2] = it_node->FastGetSolutionStepValue(REACTION_Y);
        RHS_expl[i_node * 5 + 3] = it_node->FastGetSolutionStepValue(REACTION_Z);
        RHS_expl[i_node * 5 + 4] = it_node->FastGetSolutionStepValue(REACTION_ENERGY);
    }
    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_impl, 1e-5);
}


/**
 * @brief Test the 3D explicit compressible Navier-Stokes element RHS
 * This test checks the residual calculation of the explicit compressible Navier-Stokes tetrahedra element
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicitRHS3D4N, FluidDynamicsApplicationFastSuite)
{
    // Create the test geometry
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 2);
    CreateCompressibleNavierStokesExplicit3D4NElement(r_model_part);

    // Define and set the nodal values
    const double R = 287.05; // J/(kg*K)
    const double T = 0.00247; // Temperature in K
    array_1d<double, 3> momentum;
    array_1d<double, 3> f_ext;
    f_ext[0] = 0.0; // x-volume force
    f_ext[1] = 0.0; // y-volume force
    f_ext[2] = 0.0; // z-volume force
    const double r = 0.0; // external pressure
    const double dt = 1.0e-1; // fake delta time to calculate the time derivatives

    for (auto &r_node : r_model_part.Nodes()){
        const double velocity = 2.9 * (1.0 - (r_node.Id() / 10.0));
        const double density = 1.772 * (1.0 - (r_node.Id() / 10.0));
        momentum[0] = velocity * density;
        momentum[1] = velocity * density;
        momentum[2] = velocity * density;
        const double total_energy = (density*R*T/(1.4-1)) + ((inner_prod(momentum,momentum))/(2*density));
        r_node.FastGetSolutionStepValue(DENSITY) = density;
        r_node.SetValue(DENSITY_TIME_DERIVATIVE, density / dt);
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.SetValue(MOMENTUM_TIME_DERIVATIVE, momentum / dt);
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        r_node.SetValue(TOTAL_ENERGY_TIME_DERIVATIVE, total_energy / dt);
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
    auto p_elem_expl = r_model_part.pGetElement(1);
    p_elem_expl->Initialize(r_process_info);
    p_elem_expl->AddExplicitContribution(r_process_info);

    // Check obtained RHS values
    // We check against the RHS obtained with the original implicit element without the intertial terms and the shock capturing
    std::vector<double> RHS_ref({-0.504166, -1.08898, -1.08944, -1.09073, -4.36048, -0.0996699, -0.483769, -0.533495, -0.536426, -3.08379, -0.100574, -0.512466, -0.483067, -0.527187, -3.00455, -0.100668, -0.500428, -0.512375, -0.476981, -2.93043});
    std::vector<double> RHS_expl(20);
    for (unsigned int i_node = 0; i_node < r_model_part.NumberOfNodes(); ++i_node) {
        const auto it_node = r_model_part.NodesBegin() + i_node;
        RHS_expl[i_node * 5    ] = it_node->FastGetSolutionStepValue(REACTION_DENSITY);
        RHS_expl[i_node * 5 + 1] = it_node->FastGetSolutionStepValue(REACTION_X);
        RHS_expl[i_node * 5 + 2] = it_node->FastGetSolutionStepValue(REACTION_Y);
        RHS_expl[i_node * 5 + 3] = it_node->FastGetSolutionStepValue(REACTION_Z);
        RHS_expl[i_node * 5 + 4] = it_node->FastGetSolutionStepValue(REACTION_ENERGY);
    }
    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_ref, 1e-4);
}

} // Namespace Testing
} // Namespace Kratos
