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
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_processes/embedded_skin_visualization_process.h"

namespace Kratos {
	namespace Testing {

        void SetUniqueTriangleModelPart(ModelPart& rNewModelPart){

            rNewModelPart.AddNodalSolutionStepVariable(DISTANCE);
            rNewModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rNewModelPart.AddNodalSolutionStepVariable(PRESSURE);

            // Generate a tetrahedron
            Properties::Pointer p_properties(new Properties(0));
            Node::Pointer p_point_1 = rNewModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            Node::Pointer p_point_2 = rNewModelPart.CreateNewNode(2, 10.0, 0.0, 0.0);
            Node::Pointer p_point_3 = rNewModelPart.CreateNewNode(3, 0.0, 10.0, 0.0);
            rNewModelPart.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

            // Set the nodal values
            const double p_pos = 1.0;
            const double p_neg = -1.0;
            const double t_pos = -2.0;
            const double t_neg = 2.0;
            array_1d<double, 3> v_pos = ZeroVector(3);
            array_1d<double, 3> v_neg = ZeroVector(3);
            array_1d<double, 3> d_pos = ZeroVector(3);
            array_1d<double, 3> d_neg = ZeroVector(3);
            v_pos[0] = 1.0;
            v_neg[0] = -1.0;
            d_pos[1] = 2.0;
            d_neg[1] = -2.0;
            const double level_set_height = 4.5;

            for (unsigned int i_node = 0; i_node < rNewModelPart.NumberOfNodes(); ++i_node){
                auto it_node = rNewModelPart.NodesBegin() + i_node;

                // Set DISTANCE values
                const double node_distance = (it_node->Y()) - level_set_height;
                it_node->FastGetSolutionStepValue(DISTANCE) = node_distance;

                // Set the nodal values
                if (node_distance > 0.0){
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_pos;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_pos;
                    it_node->SetValue(TEMPERATURE, t_pos);
                    it_node->SetValue(DISPLACEMENT, d_pos);
                } else {
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_neg;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_neg;
                    it_node->SetValue(TEMPERATURE, t_neg);
                    it_node->SetValue(DISPLACEMENT, d_neg);
                }
            }

            // Copy the nodal distance to ELEMENTAL_DISTANCES variable
            // (needed in case the Ausas functions are tested)
            for (unsigned int i_elem = 0; i_elem < rNewModelPart.NumberOfElements(); ++i_elem){
                auto it_elem = rNewModelPart.ElementsBegin() + i_elem;
                auto &r_geom = it_elem->GetGeometry();
                const unsigned int elem_n_nodes = r_geom.PointsNumber();
                Vector elem_dist(elem_n_nodes);
                for (unsigned int i_node = 0; i_node < elem_n_nodes; ++i_node)
                {
                    elem_dist[i_node] = r_geom[i_node].GetSolutionStepValue(DISTANCE);
                }
                it_elem->SetValue(ELEMENTAL_DISTANCES, elem_dist);
            }
        }

        void SetUniqueTetrahedronModelPart(ModelPart& rNewModelPart){

            rNewModelPart.AddNodalSolutionStepVariable(DISTANCE);
            rNewModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rNewModelPart.AddNodalSolutionStepVariable(PRESSURE);

            // Generate a tetrahedron
            Properties::Pointer p_properties(new Properties(0));
            Node::Pointer p_point_1 = rNewModelPart.CreateNewNode(1,  0.0,  0.0,  0.0);
            Node::Pointer p_point_2 = rNewModelPart.CreateNewNode(2, 10.0,  0.0,  0.0);
            Node::Pointer p_point_3 = rNewModelPart.CreateNewNode(3,  0.0, 10.0,  0.0);
            Node::Pointer p_point_4 = rNewModelPart.CreateNewNode(4,  0.0,  0.0, 10.0);
            rNewModelPart.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);

            // Set the nodal values
            const double p_pos = 1.0;
            const double p_neg = -1.0;
            const double t_pos = -2.0;
            const double t_neg = 2.0;
            array_1d<double, 3> v_pos = ZeroVector(3);
            array_1d<double, 3> v_neg = ZeroVector(3);
            array_1d<double, 3> d_pos = ZeroVector(3);
            array_1d<double, 3> d_neg = ZeroVector(3);
            v_pos[0] = 1.0;
            v_neg[0] = -1.0;
            d_pos[1] = 2.0;
            d_neg[1] = -2.0;
            const double level_set_height = 4.5;

            for (unsigned int i_node = 0; i_node < rNewModelPart.NumberOfNodes(); ++i_node){
                auto it_node = rNewModelPart.NodesBegin() + i_node;

                // Set DISTANCE values
                const double node_distance = (it_node->Z()) - level_set_height;
                it_node->FastGetSolutionStepValue(DISTANCE) = node_distance;

                // Set the nodal values
                if (node_distance > 0.0){
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_pos;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_pos;
                    it_node->SetValue(TEMPERATURE, t_pos);
                    it_node->SetValue(DISPLACEMENT, d_pos);
                } else {
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_neg;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_neg;
                    it_node->SetValue(TEMPERATURE, t_neg);
                    it_node->SetValue(DISPLACEMENT, d_neg);
                }
            }

            // Copy the nodal distance to ELEMENTAL_DISTANCES variable
            // (needed in case the Ausas functions are tested)
            for (unsigned int i_elem = 0; i_elem < rNewModelPart.NumberOfElements(); ++i_elem){
                auto it_elem = rNewModelPart.ElementsBegin() + i_elem;
                auto &r_geom = it_elem->GetGeometry();
                const unsigned int elem_n_nodes = r_geom.PointsNumber();
                Vector elem_dist(elem_n_nodes);
                for (unsigned int i_node = 0; i_node < elem_n_nodes; ++i_node){
                    elem_dist[i_node] = r_geom[i_node].GetSolutionStepValue(DISTANCE);
                }
                it_elem->SetValue(ELEMENTAL_DISTANCES, elem_dist);
            }
        }

	    /**
	     * Checks the embedded skin visualization process for a unique triangle with the standard shape functions
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessUniqueTriangleStandard, FluidDynamicsApplicationFastSuite)
		{
            // Set the test model part
            Model model;
            ModelPart& main_model_part = model.CreateModelPart("MainModelPart");
            SetUniqueTriangleModelPart(main_model_part);

            // Create the visualization model part
            ModelPart&  visualization_model_part = model.CreateModelPart("VisualizationModelPart");
            visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
            visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
            visualization_model_part.AddNodalSolutionStepVariable(PRESSURE);

            // Set the embedded skin visualization process
            Parameters visualization_settings(R"(
            {
                "level_set_type"          : "continuous",
                "shape_functions"         : "standard",
                "visualization_variables" : ["VELOCITY","PRESSURE"],
                "visualization_nonhistorical_variables" : ["TEMPERATURE","DISPLACEMENT"]
            })");

            EmbeddedSkinVisualizationProcess skin_visualization_process(
                main_model_part,
                visualization_model_part,
                visualization_settings);

            skin_visualization_process.ExecuteInitialize();
            skin_visualization_process.ExecuteBeforeSolutionLoop();
            skin_visualization_process.ExecuteInitializeSolutionStep();
            skin_visualization_process.ExecuteBeforeOutputStep();
            skin_visualization_process.ExecuteFinalizeSolutionStep();

            // Check values
            const double tolerance = 1.0e-8;
            const std::array<double, 8> expected_p{{-1.0,-1.0,1.0,-0.1,-0.1,-0.1,-0.1,-0.1}};
            const std::array<double, 8> expected_v_x{{-1.0,-1.0,1.0,-0.1,-0.1,-0.1,-0.1,-0.1}};
            const std::array<double, 8> expected_t{{2.0,2.0,-2.0,0.2,0.2,0.2,0.2,0.2}};
            const std::array<double, 8> expected_d_y{{-2.0,-2.0,2.0,-0.2,-0.2,-0.2,-0.2,-0.2}};
            unsigned int i = 0;
            for (auto& r_node : visualization_model_part.Nodes()) {
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(PRESSURE), expected_p[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_X), expected_v_x[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), expected_t[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(DISPLACEMENT_Y), expected_d_y[i], tolerance);
                i++;
            }
	    }

	    /**
	     * Checks the embedded skin visualization process for a unique triangle with the Ausas shape functions
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessUniqueTriangleAusas, FluidDynamicsApplicationFastSuite)
		{
            // Set the test model part
            Model model;
            ModelPart& main_model_part = model.CreateModelPart("MainModelPart");
            SetUniqueTriangleModelPart(main_model_part);

            // Create the visualization model part
            ModelPart& visualization_model_part = model.CreateModelPart("VisualizationModelPart");
            visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
            visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
            visualization_model_part.AddNodalSolutionStepVariable(PRESSURE);

            // Set the embedded skin visualization process
            Parameters visualization_settings(R"(
            {
                "level_set_type"          : "discontinuous",
                "shape_functions"         : "ausas",
                "visualization_variables" : ["VELOCITY","PRESSURE"],
                "visualization_nonhistorical_variables" : ["TEMPERATURE","DISPLACEMENT"]
            })");

            EmbeddedSkinVisualizationProcess skin_visualization_process(
                main_model_part,
                visualization_model_part,
                visualization_settings);

            skin_visualization_process.ExecuteInitialize();
            skin_visualization_process.ExecuteBeforeSolutionLoop();
            skin_visualization_process.ExecuteInitializeSolutionStep();
            skin_visualization_process.ExecuteBeforeOutputStep();
            skin_visualization_process.ExecuteFinalizeSolutionStep();

            // Check values
            const double tolerance = 1.0e-8;
            const std::array<double, 8> expected_p{{-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0}};
            const std::array<double, 8> expected_v_x{{-1.0,-1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0}};
            const std::array<double, 8> expected_t{{2.0,2.0,-2.0,-2.0,-2.0,2.0,2.0,2.0}};
            const std::array<double, 8> expected_d_y{{-2.0,-2.0,2.0,2.0,2.0,-2.0,-2.0,-2.0}};
            unsigned int i = 0;
            for (auto& r_node : visualization_model_part.Nodes()) {
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(PRESSURE), expected_p[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_X), expected_v_x[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), expected_t[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(DISPLACEMENT_Y), expected_d_y[i], tolerance);
                i++;
            }
	    }

        /**
	     * Checks the embedded skin visualization process for a unique triangle with the Ausas incised shape functions
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessUniqueTriangleAusasIncised, FluidDynamicsApplicationFastSuite)
		{
            // Set the test model part
            Model model;
            ModelPart& main_model_part = model.CreateModelPart("MainModelPart");
            SetUniqueTriangleModelPart(main_model_part);

            // Set ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED to test Ausas incised element
            array_1d<double,3> elem_edge_dist_extra(3,-1.0);
            elem_edge_dist_extra[1] = 0.55;
            for (unsigned int i_elem = 0; i_elem < main_model_part.NumberOfElements(); ++i_elem){
                auto it_elem = main_model_part.ElementsBegin() + i_elem;
                it_elem->SetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED, elem_edge_dist_extra);
            }

            // Create the visualization model part
            ModelPart& visualization_model_part = model.CreateModelPart("VisualizationModelPart");
            visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
            visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
            visualization_model_part.AddNodalSolutionStepVariable(PRESSURE);

            // Set the embedded skin visualization process
            Parameters visualization_settings(R"(
            {
                "level_set_type"          : "discontinuous",
                "shape_functions"         : "ausas",
                "visualization_variables" : ["VELOCITY","PRESSURE"],
                "visualization_nonhistorical_variables" : ["TEMPERATURE","DISPLACEMENT"]
            })");

            EmbeddedSkinVisualizationProcess skin_visualization_process(
                main_model_part,
                visualization_model_part,
                visualization_settings);

            skin_visualization_process.ExecuteInitialize();
            skin_visualization_process.ExecuteBeforeSolutionLoop();
            skin_visualization_process.ExecuteInitializeSolutionStep();
            skin_visualization_process.ExecuteBeforeOutputStep();
            skin_visualization_process.ExecuteFinalizeSolutionStep();

            // Check values
            const double tolerance = 1.0e-8;
            const std::array<double, 8> expected_p{{-1.0,-1.0,1.0,-0.1,1.0,-0.1,-1.0,-0.1}};
            const std::array<double, 8> expected_v_x{{-1.0,-1.0,1.0,-0.1,1.0,-0.1,-1.0,-0.1}};
            const std::array<double, 8> expected_t{{2.0,2.0,-2.0,0.2,-2.0,0.2,2.0,0.2}};
            const std::array<double, 8> expected_d_y{{-2.0,-2.0,2.0,-0.2,2.0,-0.2,-2.0,-0.2}};
            unsigned int i = 0;
            for (auto& r_node : visualization_model_part.Nodes()) {
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(PRESSURE), expected_p[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_X), expected_v_x[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), expected_t[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(DISPLACEMENT_Y), expected_d_y[i], tolerance);
                i++;
            }
	    }

	    /**
	     * Checks the embedded skin visualization process for a unique tetrahedron using the standard shape functions
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessUniqueTetrahedronStandard, FluidDynamicsApplicationFastSuite)
		{
            // Set the test model part
            Model model;
            ModelPart& main_model_part = model.CreateModelPart("MainModelPart");
            SetUniqueTetrahedronModelPart(main_model_part);

            // Create the visualization model part
            ModelPart& visualization_model_part = model.CreateModelPart("VisualizationModelPart");
            visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
            visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
            visualization_model_part.AddNodalSolutionStepVariable(PRESSURE);

            // Set the embedded skin visualization process
            Parameters visualization_settings(R"(
            {
                "level_set_type"          : "continuous",
                "shape_functions"         : "standard",
                "visualization_variables" : ["VELOCITY","PRESSURE"],
                "visualization_nonhistorical_variables" : ["TEMPERATURE","DISPLACEMENT"]
            })");

            EmbeddedSkinVisualizationProcess skin_visualization_process(
                main_model_part,
                visualization_model_part,
                visualization_settings);

            skin_visualization_process.ExecuteInitialize();
            skin_visualization_process.ExecuteBeforeSolutionLoop();
            skin_visualization_process.ExecuteInitializeSolutionStep();
            skin_visualization_process.ExecuteBeforeOutputStep();
            skin_visualization_process.ExecuteFinalizeSolutionStep();

            // Check values
            const double tolerance = 1.0e-8;
            const std::array<double, 13> expected_p{{-1.0,-1.0,-1.0,1.0,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1}};
            const std::array<double, 13> expected_v_x{{-1.0,-1.0,-1.0,1.0,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1}};
            const std::array<double, 13> expected_t{{2.0,2.0,2.0,-2.0,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2}};
            const std::array<double, 13> expected_d_y{{-2.0,-2.0,-2.0,2.0,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2}};
            unsigned int i = 0;
            for (auto& r_node : visualization_model_part.Nodes()) {
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(PRESSURE), expected_p[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_X), expected_v_x[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), expected_t[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(DISPLACEMENT_Y), expected_d_y[i], tolerance);
                i++;
            }
	    }

	    /**
	     * Checks the embedded skin visualization process for a unique tetrahedron using the Ausas shape functions
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessUniqueTetrahedronAusas, FluidDynamicsApplicationFastSuite)
		{
            // Set the test model part
            Model model;
            ModelPart& main_model_part = model.CreateModelPart ("MainModelPart");
            SetUniqueTetrahedronModelPart(main_model_part);

            // Create the visualization model part
            ModelPart& visualization_model_part = model.CreateModelPart("VisualizationModelPart");
            visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
            visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
            visualization_model_part.AddNodalSolutionStepVariable(PRESSURE);

            // Set the embedded skin visualization process
            Parameters visualization_settings(R"(
            {
                "level_set_type"          : "discontinuous",
                "shape_functions"         : "ausas",
                "visualization_variables" : ["VELOCITY","PRESSURE"],
                "visualization_nonhistorical_variables" : ["TEMPERATURE","DISPLACEMENT"]
            })");

            EmbeddedSkinVisualizationProcess skin_visualization_process(
                main_model_part,
                visualization_model_part,
                visualization_settings);

            skin_visualization_process.ExecuteInitialize();
            skin_visualization_process.ExecuteBeforeSolutionLoop();
            skin_visualization_process.ExecuteInitializeSolutionStep();
            skin_visualization_process.ExecuteBeforeOutputStep();
            skin_visualization_process.ExecuteFinalizeSolutionStep();

            // Check values
            const double tolerance = 1.0e-8;
            const std::array<double, 13> expected_p{{-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}};
            const std::array<double, 13> expected_v_x{{-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}};
            const std::array<double, 13> expected_t{{2.0,2.0,2.0,-2.0,-2.0,-2.0,-2.0,2.0,2.0,2.0,2.0,2.0,2.0}};
            const std::array<double, 13> expected_d_y{{-2.0,-2.0,-2.0,2.0,2.0,2.0,2.0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0}};
            unsigned int i = 0;
            for (auto& r_node : visualization_model_part.Nodes()) {
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(PRESSURE), expected_p[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_X), expected_v_x[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), expected_t[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(DISPLACEMENT_Y), expected_d_y[i], tolerance);
                i++;
            }
	    }

        /**
	     * Checks the embedded skin visualization process for a unique tetrahedron using the Ausas incised shape functions
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessUniqueTetrahedronAusasIncised, FluidDynamicsApplicationFastSuite)
		{
            // Set the test model part
            Model model;
            ModelPart& main_model_part = model.CreateModelPart ("MainModelPart");
            SetUniqueTetrahedronModelPart(main_model_part);

            // Set ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED to test Ausas incised element
            array_1d<double,6> elem_edge_dist_extra(6,-1.0);
            elem_edge_dist_extra[3] = 0.45;
            for (unsigned int i_elem = 0; i_elem < main_model_part.NumberOfElements(); ++i_elem){
                auto it_elem = main_model_part.ElementsBegin() + i_elem;
                it_elem->SetValue(ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED, elem_edge_dist_extra);
            }

            // Create the visualization model part
            ModelPart& visualization_model_part = model.CreateModelPart("VisualizationModelPart");
            visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
            visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
            visualization_model_part.AddNodalSolutionStepVariable(PRESSURE);

            // Set the embedded skin visualization process
            Parameters visualization_settings(R"(
            {
                "level_set_type"          : "discontinuous",
                "shape_functions"         : "ausas",
                "visualization_variables" : ["VELOCITY","PRESSURE"],
                "visualization_nonhistorical_variables" : ["TEMPERATURE","DISPLACEMENT"]
            })");

            EmbeddedSkinVisualizationProcess skin_visualization_process(
                main_model_part,
                visualization_model_part,
                visualization_settings);

            skin_visualization_process.ExecuteInitialize();
            skin_visualization_process.ExecuteBeforeSolutionLoop();
            skin_visualization_process.ExecuteInitializeSolutionStep();
            skin_visualization_process.ExecuteBeforeOutputStep();
            skin_visualization_process.ExecuteFinalizeSolutionStep();

            // Check values
            const double tolerance = 1.0e-8;
            const std::array<double, 13> expected_p{{-1.0,-1.0,-1.0,1.0,-0.1,1.0,1.0,-1.0,-0.1,-0.1,-1.0,-0.1,-1.0}};
            const std::array<double, 13> expected_v_x{{-1.0,-1.0,-1.0,1.0,-0.1,1.0,1.0,-1.0,-0.1,-0.1,-1.0,-0.1,-1.0}};
            const std::array<double, 13> expected_t{{2.0,2.0,2.0,-2.0,0.2,-2.0,-2.0,2.0,0.2,0.2,2.0,0.2,2.0}};
            const std::array<double, 13> expected_d_y{{-2.0,-2.0,-2.0,2.0,-0.2,2.0,2.0,-2.0,-0.2,-0.2,-2.0,-0.2,-2.0}};
            unsigned int i = 0;
            for (auto& r_node : visualization_model_part.Nodes()) {
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(PRESSURE), expected_p[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.FastGetSolutionStepValue(VELOCITY_X), expected_v_x[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(TEMPERATURE), expected_t[i], tolerance);
                KRATOS_EXPECT_NEAR(r_node.GetValue(DISPLACEMENT_Y), expected_d_y[i], tolerance);
                i++;
            }
	    }
    } // namespace Testing
}  // namespace Kratos.
