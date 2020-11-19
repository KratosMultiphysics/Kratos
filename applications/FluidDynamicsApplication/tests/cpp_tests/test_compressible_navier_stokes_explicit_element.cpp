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
    rModelPart.AddNodalSolutionStepVariable(MASS_SOURCE);
    rModelPart.AddNodalSolutionStepVariable(HEAT_SOURCE);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(REACTION_DENSITY);
    rModelPart.AddNodalSolutionStepVariable(MOMENTUM);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(REACTION_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    rModelPart.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    rModelPart.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);
    rModelPart.AddNodalSolutionStepVariable(CONDUCTIVITY);

    // Set the element properties
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(CONDUCTIVITY,  0.024);
    p_elem_prop->SetValue(SPECIFIC_HEAT,  722);
    p_elem_prop->SetValue(HEAT_CAPACITY_RATIO, 1.4);
    p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 2.0e-05);

    // Set process info values
    auto& r_process_info = rModelPart.GetProcessInfo();
    r_process_info[OSS_SWITCH] = false;
    r_process_info[DELTA_TIME] = 1.0e-1;
    r_process_info[TIME_INTEGRATION_THETA] = 1.0;
    r_process_info[SHOCK_CAPTURING_SWITCH] = true;

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
    rModelPart.AddNodalSolutionStepVariable(MASS_SOURCE);
    rModelPart.AddNodalSolutionStepVariable(HEAT_SOURCE);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(REACTION_DENSITY);
    rModelPart.AddNodalSolutionStepVariable(MOMENTUM);
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(REACTION_ENERGY);
    rModelPart.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
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

    // Set process info values
    auto& r_process_info = rModelPart.GetProcessInfo();
    r_process_info[OSS_SWITCH] = false;
    r_process_info[DELTA_TIME] = 1.0e-1;
    r_process_info[TIME_INTEGRATION_THETA] = 1.0;
    r_process_info[SHOCK_CAPTURING_SWITCH] = true;

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
    array_1d<double, 3> f_ext;
    array_1d<double, 3> momentum;
    f_ext[0] = 0.0; // x-volume force
    f_ext[1] = 0.0; // y-volume force
    f_ext[2] = 0.0; // z-volume force
    const double r = 0.0; // heat source
    const double mass = 0.0; // mass source

    for (auto &r_node : r_model_part.Nodes()){
        const double velocity = 2.9 * (1.0 - (r_node.Id() / 10.0));
        const double density = 1.772 * (1.0 - (r_node.Id() / 10.0));
        momentum[0] = velocity * density;
        momentum[1] = velocity * density;
        momentum[2] = 0.0;
        const double total_energy = density*R*T/(1.4-1) + inner_prod(momentum,momentum)/(2*density);
        // Set DOF values
        r_node.FastGetSolutionStepValue(DENSITY) = density;
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        // Set source terms
        r_node.FastGetSolutionStepValue(MASS_SOURCE) = mass;
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(HEAT_SOURCE) = r;
        // Set shock capturing values
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, r_node.Id() * 1e-3);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, r_node.Id() * 2e-3);
        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, r_node.Id() * 3e-3);
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
    auto p_elem_expl = r_model_part.pGetElement(1);
    p_elem_expl->Initialize(r_process_info);
    p_elem_expl->AddExplicitContribution(r_process_info);

    // Check obtained RHS values
    // We check against the RHS obtained with the original implicit element without the intertial terms and the shock capturing
    std::vector<double> RHS_ref({1.88982, 4.72962, 4.73201, 15.1248, -0.316645, -0.291717, -0.022091, 0.731414, -0.31417, -0.0713855, -0.286044, 0.618974});
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
    const double mass = 0.0; // mass source

    for (auto &r_node : r_model_part.Nodes()){
        const double velocity = 2.9 * (1.0 - (r_node.Id() / 10.0));
        const double density = 1.772 * (1.0 - (r_node.Id() / 10.0));
        momentum[0] = velocity * density;
        momentum[1] = velocity * density;
        momentum[2] = velocity * density;
        const double total_energy = (density*R*T/(1.4-1)) + ((inner_prod(momentum,momentum))/(2*density));
        // Set DOF value
        r_node.FastGetSolutionStepValue(DENSITY) = density;
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        // Set source terms
        r_node.FastGetSolutionStepValue(MASS_SOURCE) = mass;
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(HEAT_SOURCE) = r;
        // Set shock capturing values
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, r_node.Id() * 1e-3);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, r_node.Id() * 2e-3);
        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, r_node.Id() * 3e-3);
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
    auto p_elem_expl = r_model_part.pGetElement(1);
    p_elem_expl->Initialize(r_process_info);
    p_elem_expl->AddExplicitContribution(r_process_info);

    // Check obtained RHS values
    // We check against the RHS obtained with the original implicit element without the intertial terms and the shock capturing
    std::vector<double> RHS_ref({0.504157, 1.08513, 1.08448, 1.08496, 4.32898, 0.0996732, 0.48471, 0.53481, 0.538012, 3.09209, 0.100577, 0.513781, 0.484746, 0.529141, 3.01525, 0.100672, 0.502014, 0.51433, 0.479201, 2.94288});
    std::vector<double> RHS_expl(20);
    for (unsigned int i_node = 0; i_node < r_model_part.NumberOfNodes(); ++i_node) {
        const auto it_node = r_model_part.NodesBegin() + i_node;
        RHS_expl[i_node * 5    ] = it_node->FastGetSolutionStepValue(REACTION_DENSITY);
        RHS_expl[i_node * 5 + 1] = it_node->FastGetSolutionStepValue(REACTION_X);
        RHS_expl[i_node * 5 + 2] = it_node->FastGetSolutionStepValue(REACTION_Y);
        RHS_expl[i_node * 5 + 3] = it_node->FastGetSolutionStepValue(REACTION_Z);
        RHS_expl[i_node * 5 + 4] = it_node->FastGetSolutionStepValue(REACTION_ENERGY);
    }
    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_ref, 1e-5);
}

} // Namespace Testing
} // Namespace Kratos
