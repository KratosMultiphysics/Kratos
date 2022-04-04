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
#include <iomanip>


namespace Kratos {
namespace Testing {
namespace CompressibleNSConservation {

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
        auto p_cond = model_part.CreateNewCondition("LineCondition2D2N", i, node_ids, p_properties);
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

std::vector<double> AssembleReactionVector(ModelPart const& rModelPart)
{
    std::vector<double> values;
    std::size_t ndofs = rModelPart.NumberOfNodes() == 0 ? 0 : rModelPart.NumberOfNodes() * rModelPart.NodesBegin()->GetDofs().size();
    values.reserve(ndofs);

    for(auto const& r_node: rModelPart.Nodes())
    {
        values.push_back(r_node.FastGetSolutionStepValue(REACTION_DENSITY));
        values.push_back(r_node.FastGetSolutionStepValue(REACTION_X));
        values.push_back(r_node.FastGetSolutionStepValue(REACTION_Y));
        values.push_back(r_node.FastGetSolutionStepValue(REACTION_ENERGY));
    }

    return values;
}

void PrintReactions(ModelPart const& rModelPart)
{
    std::cout << "\nNode #    DENSITY           MOMENTUM_X       MOMENTUM_Y    TOTAL_ENERGY\n";

    for(auto const& r_node: rModelPart.Nodes())
    {
        std::cout << r_node.Id() << "    "
                  << std::setfill(' ') << std::right << std::setw(14)
                  << r_node.FastGetSolutionStepValue(REACTION_DENSITY) << "   "
                  << std::setfill(' ') << std::right << std::setw(14)
                  << r_node.FastGetSolutionStepValue(REACTION_X) << "   "
                  << std::setfill(' ') << std::right << std::setw(14)
                  << r_node.FastGetSolutionStepValue(REACTION_Y) << "   "
                  << std::setfill(' ') << std::right << std::setw(14)
                  << r_node.FastGetSolutionStepValue(REACTION_ENERGY) << "   \n";
    }
    std::cout << std::endl;
}
}

/**
 * @brief Test the 2D explicit compressible Navier-Stokes element and condition RHS
 * This is a conservation test.
 *
 * RIGID BODY TRANSLATION
 * N-S equations say that with:
 *  - ∇U = 0
 * then the time derivatives should be zero.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicit2D_ConservationRigidTranslation, FluidDynamicsApplicationFastSuite)
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
    const std::vector<double> RHS_ref(12, 0.0);
    std::vector<double> RHS_expl = CompressibleNSConservation::AssembleReactionVector(r_model_part);


    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_ref, 1e-4);
}

/**
 * @brief Test the 2D explicit compressible Navier-Stokes element and condition RHS
 * This is a conservation test.
 * 
 * N-S equations say that given the following conditions:
 *  - V = 0
 *  - ∇ e_total = 0
 *  - λ = 0
 * the time derivatives should be zero
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicit2D_ConservationStatic, FluidDynamicsApplicationFastSuite)
{
    // Create the test geometry
    Model model;
    ModelPart& r_model_part = GenerateModel(model);

    // Define and set the nodal values
    for (auto &r_node : r_model_part.Nodes())
    {
        const double rho = 1.2 + 0.2*r_node.X() + 0.1*r_node.Y();
        const double etot = 15;

        // Set DOF values
        r_node.FastGetSolutionStepValue(DENSITY) = rho;
        r_node.FastGetSolutionStepValue(MOMENTUM) = ZeroVector(3);
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = etot;

        r_node.FastGetSolutionStepValue(DENSITY, 1) = rho;
        r_node.FastGetSolutionStepValue(MOMENTUM, 1) = ZeroVector(3);
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1) = etot;

        // Set shock capturing values
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, 1.0);
        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 1.0);
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
    const std::vector<double> RHS_ref(12, 0.0);
    std::vector<double> RHS_expl = CompressibleNSConservation::AssembleReactionVector(r_model_part);

    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_ref, 1e-4);
}


/**
 * @brief Test the 2D explicit compressible Navier-Stokes element and condition RHS
 * This is a conservation test.
 * 
 * RIGID BODY ROTATION
 * N-S equations say that given the following conditions:
 *  - ∇ρ = 0
 *  - V = (-y, x)^T
 *  - λ = 0
 *  - p = ½ρω²r²+p0
 * the time derivatives should be zero.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressibleNavierStokesExplicit2D_ConservationRigidRotation, FluidDynamicsApplicationFastSuite)
{
    // Create the test geometry
    Model model;
    ModelPart& r_model_part = GenerateModel(model);
    
    constexpr double rho = 1.2;
    constexpr double c_v = 722.14;
    constexpr double gamma = 1.4;
    constexpr double R = c_v * (gamma - 1);
    constexpr double omega = 3;   // Angular frequency, radians per second
    constexpr double p0 = 101325; // Pressure at center of rotation, Pa

    // Define and set the nodal values
    for (auto &r_node : r_model_part.Nodes())
    {
        const double Vx = - r_node.Y();
        const double Vy = r_node.X();

        const double r2 = inner_prod(r_node, r_node);
        const double pressure = 0.5*rho*omega*omega*r2 + p0;

        const double T = pressure / rho * R;
        const double V2 = Vx*Vx + Vy*Vy;
        const double etot = rho * (0.5*V2 + c_v*T);

        // Set DOF values
        r_node.FastGetSolutionStepValue(DENSITY) = rho;
        r_node.FastGetSolutionStepValue(MOMENTUM_X) = rho * Vx;
        r_node.FastGetSolutionStepValue(MOMENTUM_Y) = rho * Vy;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = etot;

        r_node.FastGetSolutionStepValue(DENSITY, 1) = rho;
        r_node.FastGetSolutionStepValue(MOMENTUM_X, 1) = rho * Vx;
        r_node.FastGetSolutionStepValue(MOMENTUM_Y, 1) = rho * Vy;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1) = etot;

        // Set shock capturing values
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, 1.0);
        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 1.0);
    }

    // Compute explicit RHS
    const auto &r_process_info = r_model_part.GetProcessInfo();

    for(auto& r_elem: r_model_part.Elements())   { r_elem.Check(r_process_info); }
    for(auto& r_cond: r_model_part.Conditions()) { r_cond.Check(r_process_info); }

    for(auto& r_elem: r_model_part.Elements())   { r_elem.Initialize(r_process_info); }
    for(auto& r_cond: r_model_part.Conditions()) { r_cond.Initialize(r_process_info); }

    for(auto& r_cond: r_model_part.Conditions()) { r_cond.AddExplicitContribution(r_process_info); }
    for(auto& r_elem: r_model_part.Elements())   { r_elem.AddExplicitContribution(r_process_info); }

    //CompressibleNSConservation::PrintReactions(r_model_part);

    // Check obtained RHS values
    const std::vector<double> RHS_ref(12, 0.0);
    const std::vector<double> RHS_expl = CompressibleNSConservation::AssembleReactionVector(r_model_part);

    KRATOS_CHECK_VECTOR_NEAR(RHS_expl, RHS_ref, 1e-4);
}


} // Namespace Testing
} // Namespace Kratos
