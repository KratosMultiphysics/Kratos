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
 * @brief Test the 2D explicit compressible Navier-Stokes element RHS
 * This test compares the 2D explicit RHS contribution to the implicit one without the inertial terms
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicitVsImplicitRHS2D3N, FluidDynamicsApplicationFastSuite)
{
    // Supersonic test for constant variables
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 3);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_model_part.AddNodalSolutionStepVariable(REACTION_DENSITY);
    r_model_part.AddNodalSolutionStepVariable(MOMENTUM);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    r_model_part.AddNodalSolutionStepVariable(REACTION_ENERGY);
    r_model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    r_model_part.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    r_model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    r_model_part.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);
    r_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);

    // Set the element properties
    Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
    p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
    p_elem_prop->SetValue(CONDUCTIVITY,  0.0257);
    p_elem_prop->SetValue(SPECIFIC_HEAT,  718);
    p_elem_prop->SetValue(HEAT_CAPACITY_RATIO, 1.4);

    // Geometry creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    auto p_elem_expl = r_model_part.CreateNewElement("CompressibleNavierStokesExplicit2D3N", 1, elem_nodes, p_elem_prop);

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
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
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
 * @brief Test the 3D explicit compressible Navier-Stokes element RHS
 * This test compares the 3D explicit RHS contribution to the implicit one without the inertial terms
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicitVsImplicitRHS3D4N, FluidDynamicsApplicationFastSuite)
{
    // Supersonic test for constant variables
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 3);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    r_model_part.AddNodalSolutionStepVariable(DENSITY);
    r_model_part.AddNodalSolutionStepVariable(REACTION_DENSITY);
    r_model_part.AddNodalSolutionStepVariable(MOMENTUM);
    r_model_part.AddNodalSolutionStepVariable(REACTION);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    r_model_part.AddNodalSolutionStepVariable(REACTION_ENERGY);
    r_model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(DYNAMIC_VISCOSITY);
    r_model_part.AddNodalSolutionStepVariable(KINEMATIC_VISCOSITY);
    r_model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    r_model_part.AddNodalSolutionStepVariable(HEAT_CAPACITY_RATIO);
    r_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);

    // Set the element properties
    Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
    p_elem_prop->SetValue(DYNAMIC_VISCOSITY, 1.846e-05);
    p_elem_prop->SetValue(KINEMATIC_VISCOSITY, 1.568e-05);
    p_elem_prop->SetValue(CONDUCTIVITY,  0.0257);
    p_elem_prop->SetValue(SPECIFIC_HEAT,  718);
    p_elem_prop->SetValue(HEAT_CAPACITY_RATIO, 1.4);

    // Geometry creation
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3, 4};
    auto p_elem_expl = r_model_part.CreateNewElement("CompressibleNavierStokesExplicit3D4N", 1, elem_nodes, p_elem_prop);

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
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;
        r_node.FastGetSolutionStepValue(BODY_FORCE) = f_ext;
        r_node.FastGetSolutionStepValue(EXTERNAL_PRESSURE) = r;
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();
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

} // Namespace Testing
} // Namespace Kratos
