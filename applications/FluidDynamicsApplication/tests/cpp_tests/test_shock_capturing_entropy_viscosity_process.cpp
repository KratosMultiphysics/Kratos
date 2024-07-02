//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/gid_io.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "processes/structured_mesh_generator_process.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_processes/shock_capturing_entropy_viscosity_process.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos {
namespace Testing {
namespace ShockCapturingEntropyViscosityTesting{

    void SetTestModelPart(ModelPart& rNewModelPart)
    {
        // Set the problem dimensions
        rNewModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        rNewModelPart.AddNodalSolutionStepVariable(DENSITY);
        rNewModelPart.AddNodalSolutionStepVariable(PRESSURE);
        rNewModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rNewModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
        rNewModelPart.AddNodalSolutionStepVariable(NUMERICAL_ENTROPY);

        rNewModelPart.SetBufferSize(2);

        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        auto p_point_1 = Kratos::make_intrusive<Node>(1, -1.0, -1.0,  0.0);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, -1.0,  1.00,  0.0);
        auto p_point_3 = Kratos::make_intrusive<Node>(3,  1.0,  1.00,  0.0);
        auto p_point_4 = Kratos::make_intrusive<Node>(4,  1.0, -1.0,  0.0);
        Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 90,
            "element_name": "Element2D3N",
            "create_skin_sub_model_part": false
        })");
        StructuredMeshGeneratorProcess(geometry, rNewModelPart, mesher_parameters).Execute();

        // Set properties
        constexpr double c_v = 722.14;
        constexpr double gamma = 1.4;
        Properties::Pointer p_prop = Kratos::make_shared<Properties>(1);
        p_prop->SetValue(SPECIFIC_HEAT, c_v);
        p_prop->SetValue(HEAT_CAPACITY_RATIO, gamma);

        rNewModelPart.AddProperties(p_prop);
        for (Element& rElement : rNewModelPart.Elements()) {
            rElement.SetProperties(p_prop);
        }

        rNewModelPart.GetProcessInfo().GetValue(DELTA_TIME) = 1.0;
    }

    /**
     * Rankine-Hugoniot shock
     */
    void SetShockFunction(ModelPart& rNewModelPart)
    {
        for (auto& r_node : rNewModelPart.Nodes()) {
            double& r_rho  = r_node.FastGetSolutionStepValue(DENSITY);
            double& r_temp = r_node.FastGetSolutionStepValue(TEMPERATURE);
            double& r_vel  = r_node.FastGetSolutionStepValue(VELOCITY_X);
            double& r_pres = r_node.FastGetSolutionStepValue(PRESSURE);

            r_rho  = r_node.X() < 0.0 ? 1.17712 : 2.07104;
            r_temp = r_node.X() < 0.0 ?     298 : 381.699;
            r_vel  = r_node.X() < 0.0 ?     500 : 284.185;
            r_pres = r_node.X() < 0.0 ?  101325 :  228345;
        }

        rNewModelPart.CloneTimeStep(rNewModelPart.GetProcessInfo().GetValue(DELTA_TIME));
    }

    void SetSmoothFunction(ModelPart& rNewModelPart)
    {
        constexpr double total_energy = 1.07599e+06;
        constexpr double c_v = 722.14;
        constexpr double gamma = 1.4;
        constexpr double R = (gamma - 1) * c_v;

        for (auto& r_node : rNewModelPart.Nodes()) {
            double& r_rho  = r_node.FastGetSolutionStepValue(DENSITY);
            double& r_temp = r_node.FastGetSolutionStepValue(TEMPERATURE);
            double& r_vel  = r_node.FastGetSolutionStepValue(VELOCITY_X);
            double& r_pres = r_node.FastGetSolutionStepValue(PRESSURE);

            r_rho  = 5.0 - r_node.X(); // Expansion -> No shock
            r_vel  = 1.0 / r_rho;      // Respecting mass conservation
            r_temp = (total_energy / r_rho - 0.5*r_vel*r_vel) / c_v; // Constant energy
            r_pres = r_rho / (R * r_temp);
        }

        rNewModelPart.CloneTimeStep(rNewModelPart.GetProcessInfo().GetValue(DELTA_TIME));
    }
} // anonymous namespace

    /**
     * Checks the shock detection process with a smooth field with no expected shocks
     */
    KRATOS_TEST_CASE_IN_SUITE(ShockCapturingEntropyViscositySmoothField, FluidDynamicsApplicationFastSuite)
    {
        // Set the test model part
        Model model;
        auto& r_model_part = model.CreateModelPart("MainModelPart");
        ShockCapturingEntropyViscosityTesting::SetTestModelPart(r_model_part);
        ShockCapturingEntropyViscosityTesting::SetSmoothFunction(r_model_part);

        // Perform the shock detection
        Parameters sc_settings(R"(
        {
            "model_part_name" : "MainModelPart",
            "calculate_nodal_area_at_each_step" : false,
            "entropy_constant"  : 1.0,
            "energy_constant"   : 0.25,
            "artificial_mass_viscosity_Prandtl"   : 0.1,
            "artificial_conductivity_Prandtl"     : 0.1
        })");

        ShockCapturingEntropyViscosityProcess sc_process(model, sc_settings);
        sc_process.Check();
        sc_process.ExecuteInitializeSolutionStep();
        sc_process.ExecuteFinalizeSolutionStep();

        // Check values
        constexpr double tolerance = 1.0e-8;
        const auto& r_test_node = r_model_part.GetNode(5100); // Arbitrary choice

        KRATOS_EXPECT_NEAR(r_test_node.GetValue(ARTIFICIAL_MASS_DIFFUSIVITY), 0.0000467293, tolerance);
        KRATOS_EXPECT_NEAR(r_test_node.GetValue(ARTIFICIAL_CONDUCTIVITY),     0.0005555568, tolerance);
        KRATOS_EXPECT_NEAR(r_test_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY),   0.0022222273, tolerance);
        KRATOS_EXPECT_NEAR(r_test_node.GetValue(ARTIFICIAL_BULK_VISCOSITY), 0.0, tolerance);
    }

    /**
     * Checks the shock detection process with the Abgrall function
     */
    KRATOS_TEST_CASE_IN_SUITE(ShockCapturingEntropyViscosityRankineHugoniot, FluidDynamicsApplicationFastSuite)
    {
        // Set the test model part
        Model model;
        auto& r_model_part = model.CreateModelPart("MainModelPart");
        ShockCapturingEntropyViscosityTesting::SetTestModelPart(r_model_part);
        ShockCapturingEntropyViscosityTesting::SetShockFunction(r_model_part);

        // Perform the shock detection
        Parameters sc_settings(R"(
        {
            "model_part_name" : "MainModelPart",
            "calculate_nodal_area_at_each_step" : false,
            "entropy_constant"  : 1.0,
            "energy_constant"   : 0.25,
            "artificial_mass_viscosity_Prandtl"   : 0.1,
            "artificial_conductivity_Prandtl"     : 0.1
        })");

        ShockCapturingEntropyViscosityProcess sc_process(model, sc_settings);
        sc_process.Check();
        sc_process.ExecuteInitializeSolutionStep();
        sc_process.ExecuteFinalizeSolutionStep();

        // Check values
        constexpr double tolerance = 1.0e-6;
        const auto& r_shock_node = r_model_part.GetNode(4156);
        const auto& r_still_node = r_model_part.GetNode(5166);

        KRATOS_EXPECT_NEAR(r_shock_node.GetValue(ARTIFICIAL_MASS_DIFFUSIVITY), 0.1055655994, tolerance);
        KRATOS_EXPECT_NEAR(r_still_node.GetValue(ARTIFICIAL_MASS_DIFFUSIVITY), 0.0, tolerance);

        KRATOS_EXPECT_NEAR(r_shock_node.GetValue(ARTIFICIAL_CONDUCTIVITY), 0.4417240025, tolerance);
        KRATOS_EXPECT_NEAR(r_still_node.GetValue(ARTIFICIAL_CONDUCTIVITY), 0.0, tolerance);

        KRATOS_EXPECT_NEAR(r_shock_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY), 1.7668960100, tolerance);
        KRATOS_EXPECT_NEAR(r_still_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY), 0.0, tolerance);

        // This process does not add dynamic viscosity
        KRATOS_EXPECT_NEAR(r_shock_node.GetValue(ARTIFICIAL_BULK_VISCOSITY), 0.0, tolerance);
        KRATOS_EXPECT_NEAR(r_still_node.GetValue(ARTIFICIAL_BULK_VISCOSITY), 0.0, tolerance);
    }

} // namespace Testing
} // namespace Kratos.
