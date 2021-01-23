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
// #include "includes/gid_io.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "processes/calculate_nodal_area_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/find_global_nodal_neighbours_process.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_processes/shock_detection_process.h"

namespace Kratos {
namespace Testing {

    void SetTestModelPart(ModelPart& rNewModelPart)
    {
        // Set the problem dimensions
        rNewModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        // Add the variable in which the Abgrall function is stored
        rNewModelPart.AddNodalSolutionStepVariable(DENSITY);

        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, -1.0, -1.0,  0.0);
        auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, -1.0,  1.00,  0.0);
        auto p_point_3 = Kratos::make_intrusive<Node<3>>(3,  1.0,  1.00,  0.0);
        auto p_point_4 = Kratos::make_intrusive<Node<3>>(4,  1.0, -1.0,  0.0);
        Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 100,
            "element_name": "Element2D3N"
        })");
        StructuredMeshGeneratorProcess(geometry, rNewModelPart, mesher_parameters).Execute();
    }

    void SetAbgrallFunction(ModelPart& rNewModelPart)
    {
        // Set the Abgrall function values
        for (auto& r_node : rNewModelPart.Nodes()) {
            double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
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
            r_rho = x_check ? f_r + std::cos(2.0 * Globals::Pi * r_node.Y()) : f_r;
        }
    }

    /**
     * Checks the shock detection process with a smooth field with no expected shocks
     */
    KRATOS_TEST_CASE_IN_SUITE(ShockDetectionSmoothField, FluidDynamicsApplicationFastSuite)
    {
        // Set the test model part
        Model model;
        auto& r_model_part = model.CreateModelPart("MainModelPart");
        SetTestModelPart(r_model_part);

        // Set a smooth field
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(DENSITY) = r_node.X() * r_node.Y();
        }

        // Perform the shock detection
        ShockDetectionProcess shock_detection(r_model_part, DENSITY, DENSITY_GRADIENT);
        shock_detection.Execute();

        // Check values
        const double tolerance = 1.0e-8;
        KRATOS_CHECK_NEAR(r_model_part.GetNode(5100).GetValue(SHOCK_SENSOR), 0.0, tolerance);

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
		// gid_io_abgrall.WriteNodalResultsNonHistorical(NODAL_AREA, r_model_part.Nodes(), 0);
		// gid_io_abgrall.WriteNodalResultsNonHistorical(SHOCK_SENSOR, r_model_part.Nodes(), 0);
		// gid_io_abgrall.WriteNodalResultsNonHistorical(DENSITY_GRADIENT, r_model_part.Nodes(), 0);
		// gid_io_abgrall.FinalizeResults();

    }

    /**
     * Checks the shock detection process with the Abgrall function
     */
    KRATOS_TEST_CASE_IN_SUITE(ShockDetectionAbgrallFunction, FluidDynamicsApplicationFastSuite)
    {
        // Set the test model part
        Model model;
        auto& r_model_part = model.CreateModelPart("MainModelPart");
        SetTestModelPart(r_model_part);
        SetAbgrallFunction(r_model_part);

        // Perform the shock detection
        ShockDetectionProcess shock_detection(r_model_part, DENSITY, DENSITY_GRADIENT);
        shock_detection.Execute();

        // Check values
        const double tolerance = 1.0e-6;
        KRATOS_CHECK_NEAR(r_model_part.GetNode(4223).GetValue(SHOCK_SENSOR), 1.0, tolerance);
        KRATOS_CHECK_NEAR(r_model_part.GetNode(4364).GetValue(SHOCK_SENSOR), 1.0, tolerance);
        KRATOS_CHECK_NEAR(r_model_part.GetNode(7131).GetValue(SHOCK_SENSOR), 1.0, tolerance);
        KRATOS_CHECK_NEAR(r_model_part.GetNode(5309).GetValue(SHOCK_SENSOR), 0.501132, tolerance);
        KRATOS_CHECK_NEAR(r_model_part.GetNode(4810).GetValue(SHOCK_SENSOR), 0.00168451, tolerance);

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
		// gid_io_abgrall.WriteNodalResultsNonHistorical(NODAL_AREA, r_model_part.Nodes(), 0);
		// gid_io_abgrall.WriteNodalResultsNonHistorical(SHOCK_SENSOR, r_model_part.Nodes(), 0);
		// gid_io_abgrall.WriteNodalResultsNonHistorical(DENSITY_GRADIENT, r_model_part.Nodes(), 0);
		// gid_io_abgrall.FinalizeResults();

    }

} // namespace Testing
} // namespace Kratos.
