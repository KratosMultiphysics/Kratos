//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/gid_io.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "processes/structured_mesh_generator_process.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_processes/shock_capturing_physics_based_process.h"

namespace Kratos {
namespace Testing {
namespace ShockCapturingPhysicsBasedTesting{

    void SetTestModelPart(ModelPart& rNewModelPart)
    {
        // Set the problem dimensions
        rNewModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        // Add the variable in which the Abgrall function is stored
        rNewModelPart.AddNodalSolutionStepVariable(DENSITY);
        rNewModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rNewModelPart.AddNodalSolutionStepVariable(TEMPERATURE);

        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        auto p_point_1 = Kratos::make_intrusive<Node>(1, -1.0, -1.0,  0.0);
        auto p_point_2 = Kratos::make_intrusive<Node>(2, -1.0,  1.0,  0.0);
        auto p_point_3 = Kratos::make_intrusive<Node>(3,  1.0,  1.0,  0.0);
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
        const double c_v = 1000;
        const double gamma = 1.4;
        Properties::Pointer p_prop = Kratos::make_shared<Properties>(1);
        p_prop->SetValue(SPECIFIC_HEAT, c_v);
        p_prop->SetValue(HEAT_CAPACITY_RATIO, gamma);

        rNewModelPart.AddProperties(p_prop);
        for (Element& rElement : rNewModelPart.Elements()) {
            rElement.SetProperties(p_prop);
        }
    }

    void SetAbgrallFunction(
        ModelPart& rNewModelPart,
        const double AvoidNegativeValues = true)
    {
        // Set the Abgrall function values
        double min_aux_val = std::numeric_limits<double>::max();
        for (auto& r_node : rNewModelPart.Nodes()) {
            double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            double& r_temp = r_node.FastGetSolutionStepValue(TEMPERATURE);
            array_1d<double,3>& r_vel = r_node.FastGetSolutionStepValue(VELOCITY);
            const bool x_check = r_node.X() > std::cos(Globals::Pi * r_node.Y()) / 2.0;
            const double y_val = r_node.Y() / std::tan(std::sqrt(Globals::Pi / 2.0));
            const double r = x_check ? r_node.X() + y_val : r_node.X() - y_val;
            double f_r;
            const double one_third = 1.0 / 3.0;
            if (std::abs(r) < one_third) {
                f_r = std::abs(std::sin(2.0 * Globals::Pi * r));
            } else if (r >= one_third) {
                f_r = 2.0 * r - 1 + (1.0 / 6.0) * std::sin(3.0 * Globals::Pi * r);
            } else {
                f_r = - r * std::sin(1.5 * Globals::Pi * std::pow(r,2));
            }
            const double aux_val = x_check ? f_r + std::cos(2.0 * Globals::Pi * r_node.Y()) : f_r;
            r_rho = aux_val;
            r_temp = aux_val;
            r_vel[0] = aux_val;
            r_vel[1] = aux_val;
            if (aux_val < min_aux_val) {
                min_aux_val = aux_val;
            }
        }

        // Avoid negative (non-physical) values in the Abgrall function
        const double epsilon = 1.0e-7;
        const double neg_threshold = std::abs(min_aux_val) + epsilon;
        if (AvoidNegativeValues) {
            for (auto& r_node : rNewModelPart.Nodes()) {
                double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
                double& r_temp = r_node.FastGetSolutionStepValue(TEMPERATURE);
                array_1d<double,3>& r_vel = r_node.FastGetSolutionStepValue(VELOCITY);
                r_rho += neg_threshold;
                r_temp += neg_threshold;
                r_vel[0] += neg_threshold;
                r_vel[1] += neg_threshold;
            }
        }
    }
}

    /**
     * Checks the shock detection process with a smooth field with no expected shocks
     */
    KRATOS_TEST_CASE_IN_SUITE(ShockCapturingPhysicsBasedSmoothField, FluidDynamicsApplicationFastSuite)
    {
        // Set the test model part
        Model model;
        auto& r_model_part = model.CreateModelPart("MainModelPart");
        ShockCapturingPhysicsBasedTesting::SetTestModelPart(r_model_part);

        // Set a smooth field
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(DENSITY) = 1.0 + std::abs(r_node.X() * r_node.Y());
            r_node.FastGetSolutionStepValue(TEMPERATURE) = 1.0 + 2.0 * std::abs(r_node.X() * r_node.Y());
            r_node.FastGetSolutionStepValue(VELOCITY_X) = 1.0 + std::pow(r_node.X(),2) + std::pow(r_node.Y(),2);
        }

        // Perform the shock detection
        Parameters sc_settings(R"(
        {
            "model_part_name" : "MainModelPart",
            "shock_sensor" : true,
            "shear_sensor" : true,
            "thermal_sensor" : true,
            "thermally_coupled_formulation" : true
        })");
        ShockCapturingPhysicsBasedProcess sc_process(model, sc_settings);
        sc_process.Execute();

        // Check values
        const double tolerance = 1.0e-8;
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(5100).GetValue(ARTIFICIAL_CONDUCTIVITY), 0.000107408, tolerance);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(5100).GetValue(ARTIFICIAL_BULK_VISCOSITY), 0.000902773, tolerance);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(5100).GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY), 6.45302e-05, tolerance);

        // GidIO<> gid_io_abgrall(
        //     "/home/rzorrilla/Desktop/abgrall_function_smooth_field",
        //     GiD_PostAscii,
        //     SingleFile,
        //     WriteDeformed,
        //     WriteConditions);
		// gid_io_abgrall.InitializeMesh(0.0);
		// gid_io_abgrall.WriteMesh(r_model_part.GetMesh());
		// gid_io_abgrall.FinalizeMesh();
		// gid_io_abgrall.InitializeResults(0, r_model_part.GetMesh());
        // gid_io_abgrall.WriteNodalResults(DENSITY, r_model_part.Nodes(), 0, 0);
        // gid_io_abgrall.WriteNodalResults(VELOCITY, r_model_part.Nodes(), 0, 0);
        // gid_io_abgrall.WriteNodalResults(TEMPERATURE, r_model_part.Nodes(), 0, 0);
        // gid_io_abgrall.WriteNodalResultsNonHistorical(ARTIFICIAL_CONDUCTIVITY, r_model_part.Nodes(), 0);
        // gid_io_abgrall.WriteNodalResultsNonHistorical(ARTIFICIAL_BULK_VISCOSITY, r_model_part.Nodes(), 0);
        // gid_io_abgrall.WriteNodalResultsNonHistorical(ARTIFICIAL_DYNAMIC_VISCOSITY, r_model_part.Nodes(), 0);
        // gid_io_abgrall.FinalizeResults();

    }

    /**
     * Checks the shock detection process with the Abgrall function
     */
    KRATOS_TEST_CASE_IN_SUITE(ShockCapturingPhysicsBasedAbgrallFunction, FluidDynamicsApplicationFastSuite)
    {
        // Set the test model part
        Model model;
        auto& r_model_part = model.CreateModelPart("MainModelPart");
        ShockCapturingPhysicsBasedTesting::SetTestModelPart(r_model_part);
        ShockCapturingPhysicsBasedTesting::SetAbgrallFunction(r_model_part);

        // Perform the shock detection
        Parameters sc_settings(R"(
        {
            "model_part_name" : "MainModelPart",
            "shock_sensor" : true,
            "shear_sensor" : true,
            "thermal_sensor" : true,
            "thermally_coupled_formulation" : true
        })");
        ShockCapturingPhysicsBasedProcess sc_process(model, sc_settings);
        sc_process.Execute();

        // Check values
        const double tolerance = 1.0e-8;
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(4223).GetValue(ARTIFICIAL_CONDUCTIVITY), 0.000188616, tolerance);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(4364).GetValue(ARTIFICIAL_CONDUCTIVITY), 0.000303507, tolerance);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(7131).GetValue(ARTIFICIAL_BULK_VISCOSITY), 0.000283975, tolerance);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(5309).GetValue(ARTIFICIAL_BULK_VISCOSITY), 0.000541887, tolerance);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(4810).GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY), 3.08545e-05, tolerance);
        KRATOS_EXPECT_NEAR(r_model_part.GetNode(6000).GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY), 0.000177703, tolerance);

        // GidIO<> gid_io_abgrall(
        //     "/home/rzorrilla/Desktop/abgrall_function_shock_detection",
        //     GiD_PostAscii,
        //     SingleFile,
        //     WriteDeformed,
        //     WriteConditions);
		// gid_io_abgrall.InitializeMesh(0.0);
		// gid_io_abgrall.WriteMesh(r_model_part.GetMesh());
		// gid_io_abgrall.FinalizeMesh();
		// gid_io_abgrall.InitializeResults(0, r_model_part.GetMesh());
		// gid_io_abgrall.WriteNodalResults(DENSITY, r_model_part.Nodes(), 0, 0);
		// gid_io_abgrall.WriteNodalResults(VELOCITY, r_model_part.Nodes(), 0, 0);
		// gid_io_abgrall.WriteNodalResults(TEMPERATURE, r_model_part.Nodes(), 0, 0);
		// gid_io_abgrall.WriteNodalResultsNonHistorical(ARTIFICIAL_CONDUCTIVITY, r_model_part.Nodes(), 0);
		// gid_io_abgrall.WriteNodalResultsNonHistorical(ARTIFICIAL_BULK_VISCOSITY, r_model_part.Nodes(), 0);
		// gid_io_abgrall.WriteNodalResultsNonHistorical(ARTIFICIAL_DYNAMIC_VISCOSITY, r_model_part.Nodes(), 0);
		// gid_io_abgrall.FinalizeResults();

    }

} // namespace Testing
} // namespace Kratos.
