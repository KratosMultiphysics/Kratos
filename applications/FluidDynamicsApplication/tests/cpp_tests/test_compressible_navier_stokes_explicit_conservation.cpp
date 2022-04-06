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
#include <iomanip>
#include <sstream>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/properties.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "testing/testing.h"
#include "input_output/logger.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos {
namespace Testing {
namespace CompressibleNSConservation {


ModelPart& GenerateModel(Model& rModel, const std::string& rConditionName)
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
    p_properties->SetValue(CONDUCTIVITY,  0.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 0.0);

    p_properties->SetValue(SPECIFIC_HEAT,  722.14);
    p_properties->SetValue(HEAT_CAPACITY_RATIO, 1.4);

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
        auto p_cond = model_part.CreateNewCondition(rConditionName, i, node_ids, p_properties);
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
 * Debugging tool that prints a table with the value of all DOFs in each node
 */
void PrintReactions(ModelPart const& rModelPart)
{
    std::stringstream ss;
    ss << "\nNode #    DENSITY           MOMENTUM_X       MOMENTUM_Y    TOTAL_ENERGY\n";

    for(auto const& r_node: rModelPart.Nodes())
    {
        ss << r_node.Id() << "    "
           << std::setfill(' ') << std::right << std::setw(14)
           << r_node.FastGetSolutionStepValue(REACTION_DENSITY) << "   "
           << std::setfill(' ') << std::right << std::setw(14)
           << r_node.FastGetSolutionStepValue(REACTION_X) << "   "
           << std::setfill(' ') << std::right << std::setw(14)
           << r_node.FastGetSolutionStepValue(REACTION_Y) << "   "
           << std::setfill(' ') << std::right << std::setw(14)
           << r_node.FastGetSolutionStepValue(REACTION_ENERGY) << "   \n";
    }
    ss << std::endl;
    Logger("") << ss.str();
}

void SetDofValues(
    ModelPart& rModelPart,
    double(*Rho)(array_1d<double, 3> const&),
    array_1d<double, 3>(*Momentum)(array_1d<double, 3> const&),
    double(*TotalEnergy)(array_1d<double, 3> const&))
{
    for (auto &r_node : rModelPart.Nodes())
    {
        const double rho = Rho(r_node);
        const array_1d<double, 3> momentum = Momentum(r_node);
        const double etot = TotalEnergy(r_node);

        // Set DOF values
        r_node.FastGetSolutionStepValue(DENSITY) = rho;
        r_node.FastGetSolutionStepValue(MOMENTUM) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY) = etot;

        r_node.FastGetSolutionStepValue(DENSITY, 1) = rho;
        r_node.FastGetSolutionStepValue(MOMENTUM, 1) = momentum;
        r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1) = etot;
    }
}

void SetViscosities(
    ModelPart& rModelPart,
    const double DynamicViscosity, 
    const double BulkViscosity, 
    const double Conductivity)
{
    for (auto &r_node : rModelPart.Nodes())
    {
        r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, DynamicViscosity);
        r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, BulkViscosity);
        r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, Conductivity);
    }
}

void SetEulerFluxes(
    ModelPart& rModelPart,
    double(*Rho)(array_1d<double, 3> const&),
    array_1d<double, 3>(*Momentum)(array_1d<double, 3> const&),
    double(*TotalEnergy)(array_1d<double, 3> const&))
{
    for (auto& r_condition: rModelPart.Conditions())
    {
        const array_1d<double, 3> X = r_condition.GetGeometry().Center();
        const array_1d<double, 3> n = r_condition.GetGeometry().UnitNormal(0);

        const double rho = Rho(X);
        const array_1d<double, 3> V = Momentum(X) / rho;
        const double etot = TotalEnergy(X);

        constexpr double gamma = 1.4;
        const double p = (gamma - 1) * (etot - rho * inner_prod(V,V)/2);
        const double Vn = inner_prod(V, n);

        r_condition.SetValue(DENSITY_FLUX, rho * Vn);
        r_condition.SetValue(MOMENTUM_FLUX, rho * V * Vn + p * n);
        r_condition.SetValue(TOTAL_ENERGY_FLUX, (etot + p)*Vn);
    }
}

std::vector<double> Assemble(ModelPart& rModelPart, const bool debug_prints = false)
{
    const auto &r_process_info = rModelPart.GetProcessInfo();

    void(*print_reactions)(const ModelPart&) = [](const ModelPart&) { };
    if(debug_prints) {
        print_reactions = [](const ModelPart& r_model_part) { PrintReactions(r_model_part); };
    }

    // Computing reactions
    for(auto& r_elem: rModelPart.Elements())   { r_elem.Check(r_process_info); }
    for(auto& r_cond: rModelPart.Conditions()) { r_cond.Check(r_process_info); }

    for(auto& r_elem: rModelPart.Elements())   { r_elem.Initialize(r_process_info); }
    for(auto& r_cond: rModelPart.Conditions()) { r_cond.Initialize(r_process_info); }

    for(auto& r_cond: rModelPart.Conditions()) { r_cond.AddExplicitContribution(r_process_info); print_reactions(rModelPart); }
    for(auto& r_elem: rModelPart.Elements())   { r_elem.AddExplicitContribution(r_process_info); print_reactions(rModelPart); }
    
    // Assembling DOF vector
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
    using namespace CompressibleNSConservation;

    // Create the test geometry
    Model model;
    ModelPart& r_model_part = GenerateModel(model, "CompressibleNavierStokesExplicitNeumannCondition2D2N");

    // Define and set the nodal values
    constexpr double gamma = 1.4;
    constexpr double p = 101325;
    constexpr double u = 5;
    constexpr double density = 1.2;

    const auto rho  = [](array_1d<double, 3> const&) { return density; };
    const auto mom  = [](array_1d<double, 3> const&) { return array_1d<double,3>{density * u, 0, 0}; };
    const auto etot = [](array_1d<double, 3> const&) { return 0.5*density*u*u + p / (gamma - 1); };

    SetDofValues(r_model_part, rho, mom, etot);
    SetViscosities(r_model_part, 1e-3, 2e-3, 3e-3);
    SetEulerFluxes(r_model_part, rho, mom, etot);
    const auto rhs = Assemble(r_model_part);

    // Check obtained RHS values
    const std::vector<double> reference(12, 0.0);
    KRATOS_CHECK_VECTOR_NEAR(rhs, reference, 1e-4);
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
    using namespace CompressibleNSConservation;

    // Create the test geometry
    Model model;
    ModelPart& r_model_part = GenerateModel(model, "CompressibleNavierStokesExplicitNeumannCondition2D2N");

    constexpr double p = 101325;
    constexpr double gamma = 1.4;

    const auto rho  = [](array_1d<double, 3> const& X) { return 1.2 + 0.2*X[0] + 0.1*X[1]; };
    const auto mom  = [](array_1d<double, 3> const&)   { return array_1d<double,3>({0, 0, 0}); };
    const auto etot = [](array_1d<double, 3> const& X) { return p / (gamma - 1); };

    SetDofValues(r_model_part, rho, mom, etot);
    SetViscosities(r_model_part, 1e-3, 2e-3, 0.0);
    SetEulerFluxes(r_model_part, rho, mom, etot);
    const auto rhs = Assemble(r_model_part);

    // Check obtained RHS values
    const std::vector<double> reference(12, 0.0);
    KRATOS_CHECK_VECTOR_NEAR(rhs, reference, 1e-4);
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
    using namespace CompressibleNSConservation;

    // Create the test geometry
    Model model;
    ModelPart& r_model_part = GenerateModel(model, "CompressibleNavierStokesExplicitNeumannCondition2D2N");
    
    constexpr double density = 1.2;
    constexpr double gamma = 1.4;
    constexpr double omega = 3;   // Angular frequency, radians per second
    constexpr double p0 = 101325; // Pressure at center of rotation, Pa

    const auto rho  = [](array_1d<double, 3> const& X) { return density; };
    const auto mom  = [](array_1d<double, 3> const& X) { return array_1d<double, 3>{- density*omega*X[1], density*omega*X[0], 0.0}; };
    const auto etot = [](array_1d<double, 3> const& X) { return 0.5*density*omega*omega*inner_prod(X,X) + p0/(gamma-1); };

    SetDofValues(r_model_part, rho, mom, etot);
    SetViscosities(r_model_part, 1e-3, 2e-3, 0.0);
    SetEulerFluxes(r_model_part, rho, mom, etot);
    const auto rhs = Assemble(r_model_part);

    // Check obtained RHS values
    const std::vector<double> reference(12, 0.0);
    KRATOS_CHECK_VECTOR_NEAR(rhs, reference, 1e-4);
}


} // Namespace Testing
} // Namespace Kratos
