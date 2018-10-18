//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
            Node<3>::Pointer p_point_1 = rNewModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            Node<3>::Pointer p_point_2 = rNewModelPart.CreateNewNode(2, 10.0, 0.0, 0.0);
            Node<3>::Pointer p_point_3 = rNewModelPart.CreateNewNode(3, 0.0, 10.0, 0.0);
            rNewModelPart.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties);

            // Set the nodal values
            const double p_pos = 1.0;
            const double p_neg = -1.0;
            array_1d<double, 3> v_pos = ZeroVector(3);
            array_1d<double, 3> v_neg = ZeroVector(3);
            v_pos[0] = 1.0;
            v_neg[0] = -1.0;
            const double level_set_height = 4.5;

            for (unsigned int i_node = 0; i_node < rNewModelPart.NumberOfNodes(); ++i_node){
                auto it_node = rNewModelPart.NodesBegin() + i_node;

                // Set DISTANCE values
                const double node_distance = (it_node->Y()) - level_set_height;
                it_node->FastGetSolutionStepValue(DISTANCE) = node_distance;

                // Set the VELOCITY and PRESSURE values
                if (node_distance > 0.0){
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_pos;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_pos;
                } else {
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_neg;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_neg;
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
            Node<3>::Pointer p_point_1 = rNewModelPart.CreateNewNode(1,  0.0,  0.0,  0.0);
            Node<3>::Pointer p_point_2 = rNewModelPart.CreateNewNode(2, 10.0,  0.0,  0.0);
            Node<3>::Pointer p_point_3 = rNewModelPart.CreateNewNode(3,  0.0, 10.0,  0.0);
            Node<3>::Pointer p_point_4 = rNewModelPart.CreateNewNode(4,  0.0,  0.0, 10.0);
            rNewModelPart.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties);

            // Set the nodal values
            const double p_pos = 1.0;
            const double p_neg = -1.0;
            array_1d<double,3> v_pos = ZeroVector(3);
            array_1d<double,3> v_neg = ZeroVector(3);
            v_pos[0] = 1.0;
            v_neg[0] = -1.0;
            const double level_set_height = 4.5;

            for (unsigned int i_node = 0; i_node < rNewModelPart.NumberOfNodes(); ++i_node){
                auto it_node = rNewModelPart.NodesBegin() + i_node;

                // Set DISTANCE values
                const double node_distance = (it_node->Z()) - level_set_height;
                it_node->FastGetSolutionStepValue(DISTANCE) = node_distance;

                // Set the VELOCITY and PRESSURE values
                if (node_distance > 0.0){
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_pos;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_pos;
                } else {
                    it_node->FastGetSolutionStepValue(VELOCITY) = v_neg;
                    it_node->FastGetSolutionStepValue(PRESSURE) = p_neg;
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
                "shape_functions"         : "standard",
                "visualization_variables" : ["VELOCITY","PRESSURE"]
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
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
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
                "shape_functions"         : "ausas",
                "visualization_variables" : ["VELOCITY","PRESSURE"]
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
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
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
                "shape_functions"         : "standard",
                "visualization_variables" : ["VELOCITY","PRESSURE"]
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
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(8).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(9).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(10).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(11).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(12).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(13).FastGetSolutionStepValue(PRESSURE), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(8).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(9).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(10).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(11).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(12).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(13).FastGetSolutionStepValue(VELOCITY_X), -0.1, 1e-8);
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
                "shape_functions"         : "ausas",
                "visualization_variables" : ["VELOCITY","PRESSURE"]
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
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(PRESSURE), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(8).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(9).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(10).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(11).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(12).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(13).FastGetSolutionStepValue(PRESSURE), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(1).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(2).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(3).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(4).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(5).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(6).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(7).FastGetSolutionStepValue(VELOCITY_X), 1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(8).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(9).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(10).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(11).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(12).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
            KRATOS_CHECK_NEAR(visualization_model_part.GetNode(13).FastGetSolutionStepValue(VELOCITY_X), -1.0, 1e-8);
	    }
    } // namespace Testing
}  // namespace Kratos.
