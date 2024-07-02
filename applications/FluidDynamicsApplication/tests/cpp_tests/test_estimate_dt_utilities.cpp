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
#include "containers/model.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/estimate_dt_utilities.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos {
namespace Testing  {

namespace Internals {

void TestEstimateDtUtilitiesInitializeModelPart(
    ModelPart& rModelPart,
    const double DeltaTime)
{

    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    Properties::Pointer p_properties = rModelPart.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, 0.1);
    p_properties->SetValue(CONDUCTIVITY, 2.0);
    p_properties->SetValue(SPECIFIC_HEAT, 1.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0);

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
    rModelPart.CreateNewNode(4, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> element_nodes_1{1, 2, 3};
    std::vector<ModelPart::IndexType> element_nodes_2{3, 2, 4};
    rModelPart.CreateNewElement("Element2D3N", 1, element_nodes_1, p_properties);
    rModelPart.CreateNewElement("Element2D3N", 2, element_nodes_2, p_properties);

    // Set a fake current delta time
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, DeltaTime);

    // Set nodal data
    for (auto& rNode : rModelPart.Nodes()) {
        rNode.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, rNode.Id());
        rNode.SetValue(ARTIFICIAL_CONDUCTIVITY, 2.0 * rNode.Id());
        rNode.SetValue(SOUND_VELOCITY, 340.0);
        rNode.FastGetSolutionStepValue(DENSITY) = rNode.Id() / 10.0;
        rNode.FastGetSolutionStepValue(VELOCITY_X) = rNode.Id() * rNode.X();
        rNode.FastGetSolutionStepValue(VELOCITY_Y) = rNode.Id() * rNode.Y();
    }
}

} // namespace internals

KRATOS_TEST_CASE_IN_SUITE(EstimateDtUtilitiesEstimateDt, FluidDynamicsApplicationFastSuite)
{
    // Set an extremely large current delta time to obtain a large CFL number
    const double current_dt = 1.0;

    // Create the test model part
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    Internals::TestEstimateDtUtilitiesInitializeModelPart(r_model_part, current_dt);

    // Estimate the delta time
    Parameters estimate_dt_settings = Parameters(R"({
        "automatic_time_step"   : true,
        "CFL_number"            : 1.0,
        "minimum_delta_time"    : 1e-4,
        "maximum_delta_time"    : 1e+1
    })");
    const auto estimate_dt_utility = EstimateDtUtility(r_model_part, estimate_dt_settings);
    const double obtained_dt = estimate_dt_utility.EstimateDt();

    // Check results
    const double tolerance = 2.0e-6;
    const double expected_dt = 0.126211;
    KRATOS_EXPECT_NEAR(expected_dt, obtained_dt, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(EstimateDtUtilitiesEstimateDtCompressibleFlow, FluidDynamicsApplicationFastSuite)
{
    // Set an extremely large current delta time to obtain a large CFL number
    const double current_dt = 1.0;

    // Create the test model part
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    Internals::TestEstimateDtUtilitiesInitializeModelPart(r_model_part, current_dt);

    // Estimate the delta time
    Parameters estimate_dt_settings = Parameters(R"({
        "automatic_time_step"             : true,
        "CFL_number"                      : 1.0,
        "Viscous_Fourier_number"          : 1.0,
        "Thermal_Fourier_number"          : 1.0,
        "minimum_delta_time"              : 1e-4,
        "maximum_delta_time"              : 1e+1,
        "consider_artificial_diffusion"   : true,
        "nodal_density_formulation"       : true,
        "consider_compressibility_in_CFL" : true
    })");
    const auto estimate_dt_utility = EstimateDtUtility(r_model_part, estimate_dt_settings);
    const double obtained_dt = estimate_dt_utility.EstimateDt();

    // Check results
    const double tolerance = 1.0e-6;
    const double expected_dt = 0.0013017675;
    KRATOS_EXPECT_NEAR(expected_dt, obtained_dt, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(EstimateDtUtilitiesEstimateDtThermal, FluidDynamicsApplicationFastSuite)
{
    // Set an extremely large current delta time to obtain a large CFL number
    const double current_dt = 1.0;

    // Create the test model part
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("TestModelPart");
    Internals::TestEstimateDtUtilitiesInitializeModelPart(r_model_part, current_dt);

    // Estimate the delta time
    Parameters estimate_dt_settings = Parameters(R"({
        "automatic_time_step"           : true,
        "CFL_number"                    : 1.0,
        "Viscous_Fourier_number"        : 0.0,
        "Thermal_Fourier_number"        : 1.0,
        "minimum_delta_time"            : 1e-4,
        "maximum_delta_time"            : 1e+1,
        "consider_artificial_diffusion" : false,
        "nodal_density_formulation"     : false
    })");
    const auto estimate_dt_utility = EstimateDtUtility(r_model_part, estimate_dt_settings);
    const double obtained_dt = estimate_dt_utility.EstimateDt();

    // Check results
    const double tolerance = 1.0e-6;
    const double expected_dt = 0.01;
    KRATOS_EXPECT_NEAR(expected_dt, obtained_dt, tolerance);
}

}
}