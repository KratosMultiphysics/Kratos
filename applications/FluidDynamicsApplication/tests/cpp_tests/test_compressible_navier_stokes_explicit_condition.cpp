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

// External includes

// Project includes
#include "containers/array_1d.h"
#include "testing/testing.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "containers/model.h"
#include "custom_elements/compressible_navier_stokes_explicit.h"
#include "custom_conditions/compressible_navier_stokes_explicit_condition.h"

namespace Kratos {
namespace Testing {

ModelPart& GenerateModel(Model& rModel)
{
    auto& model_part = rModel.CreateModelPart("main", 2);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(MASS_SOURCE);
    model_part.AddNodalSolutionStepVariable(HEAT_SOURCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(REACTION_DENSITY);
    model_part.AddNodalSolutionStepVariable(MOMENTUM);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(TOTAL_ENERGY);
    model_part.AddNodalSolutionStepVariable(REACTION_ENERGY);

    // Set the element properties
    Properties::Pointer p_properties = model_part.CreateNewProperties(0);
    p_properties->SetValue(CONDUCTIVITY,  0.024);
    p_properties->SetValue(SPECIFIC_HEAT,  722.14);
    p_properties->SetValue(HEAT_CAPACITY_RATIO, 1.4);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 2.0e-05);

    // Set process info values
    auto& r_process_info = model_part.GetProcessInfo();
    r_process_info[OSS_SWITCH] = false;
    r_process_info[DELTA_TIME] = 1.0e-1;
    r_process_info[TIME_INTEGRATION_THETA] = 1.0;
    r_process_info[SHOCK_CAPTURING_SWITCH] = true;

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3};
    auto p_elem = model_part.CreateNewElement("CompressibleNavierStokesExplicit2D3N", 1, elem_nodes, p_properties);

    auto neighbour_elems = NEIGHBOUR_ELEMENTS.Zero();
    neighbour_elems.push_back(p_elem);

    for(std::size_t i=1; i<=3; ++i)
    {
        std::vector<std::size_t> node_ids = {i, i%3 + 1};
        auto p_cond = model_part.CreateNewCondition("CompressibleNavierStokesExplicitCondition2D2N", i, node_ids, p_properties);
        p_cond->SetValue(NEIGHBOUR_ELEMENTS, neighbour_elems);
    }

    for(auto& r_node: model_part.Nodes())
    {
        r_node.AddDof(DENSITY);
        r_node.AddDof(MOMENTUM_X);
        r_node.AddDof(MOMENTUM_Y);
        r_node.AddDof(TOTAL_ENERGY);
    }

    return model_part;
}

/**
 * @brief Test the 2D explicit compressible Navier-Stokes element and condition RHS
 * This is a conservation test.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicitConditionRHS2D2N, FluidDynamicsApplicationFastSuite)
{
    // Create the test geometry
    Model model;
    ModelPart& r_model_part = GenerateModel(model);

    // Define and set the nodal values
    constexpr double cv = 722.14;
    constexpr double T = 300; // Temperature in K
    
    constexpr double density = 1.2;
    constexpr double velocity = 5.0;
    constexpr double total_energy = density * (0.5*velocity*velocity + cv*T);

    const array_1d<double, 3> momentum {density * velocity, 0.0, 0.0};

    for (auto &r_node : r_model_part.Nodes())
    {
        // Set DOF values
        r_node.FastGetSolutionStepValue(DENSITY) = density;
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = total_energy;

        r_node.FastGetSolutionStepValue(DENSITY, 1) = density;
        r_node.FastGetSolutionStepValue(MOMENTUM, 1) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1) = total_energy;

        // Set shock capturing values
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, r_node.Id() * 1e-3);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, r_node.Id() * 2e-3);
        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, r_node.Id() * 3e-3);
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();

    for(auto& r_elem: r_model_part.Elements())   { r_elem.Check(r_process_info); }
    for(auto& r_cond: r_model_part.Conditions()) { r_cond.Check(r_process_info); }

    for(auto& r_elem: r_model_part.Elements())   { r_elem.Initialize(r_process_info); }
    for(auto& r_cond: r_model_part.Conditions()) { r_cond.Initialize(r_process_info); }

    for(auto& r_elem: r_model_part.Elements())   { r_elem.AddExplicitContribution(r_process_info); }
    for(auto& r_cond: r_model_part.Conditions()) { r_cond.AddExplicitContribution(r_process_info); }

    // Check obtained RHS values
    std::vector<double> RHS_ref(12, 0.0);
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

} // Namespace Testing
} // Namespace Kratos
