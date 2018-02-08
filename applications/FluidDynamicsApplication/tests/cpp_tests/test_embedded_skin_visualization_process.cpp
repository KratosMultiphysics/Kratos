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
#include "includes/gid_io.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_processes/embedded_skin_visualization_process.h"

namespace Kratos {
	namespace Testing {

	    /** 
	     * Checks the embedded skin visualization process for triangular meshes.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessTriangle, FluidDynamicsApplicationFastSuite)
		{

            // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
            Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 10.00,  0.00, 0.00);
            Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.00,   0.00, 0.00);
            Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 0.00,  10.00, 0.00);
            Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 10.00, 10.00, 0.00);

            Quadrilateral2D4<Node<3>> square_domain(p_point_1, p_point_2, p_point_3, p_point_4);

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions" :  5,
                "element_name"        : "Element2D3N"
            })");

            ModelPart main_model_part("MainModelPart");
            main_model_part.AddNodalSolutionStepVariable(DISTANCE);
            main_model_part.AddNodalSolutionStepVariable(VELOCITY);
            main_model_part.AddNodalSolutionStepVariable(PRESSURE);
            StructuredMeshGeneratorProcess(square_domain, main_model_part, mesher_parameters).Execute();

            // Set the nodal values
            const double level_set_height = 4.5;
            for (unsigned int i_node = 0; i_node < main_model_part.NumberOfNodes(); ++i_node){
                auto it_node = main_model_part.NodesBegin() + i_node;

                // Set DISTANCE values
                const double node_distance = (it_node->Y()) - level_set_height;
                it_node->FastGetSolutionStepValue(DISTANCE) = node_distance;

                // Set the VELOCITY and PRESSURE values
                double p_node;
                array_1d<double, 3> v_node = ZeroVector(3);
                if (node_distance > 0.0){
                    v_node[0] = it_node->Y();
                    v_node[1] = std::pow(it_node->Y(), 2);
                    p_node = (it_node->X())*(it_node->Y());
                } else {
                    v_node[0] = -(it_node->Y());
                    v_node[1] = -(std::pow(it_node->Y(), 2));
                    p_node = -(it_node->X())*(it_node->Y());
                }

                it_node->FastGetSolutionStepValue(VELOCITY) = v_node;
                it_node->FastGetSolutionStepValue(PRESSURE) = p_node;

            }

            // Create the visualization model part
            ModelPart visualization_model_part("VisualizationModelPart");
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

            // GidIO<> gid_io_fluid("/home/rzorrilla/Kratos/tests/visualizaton_model_part_2d", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io_fluid.InitializeMesh(0);
            // gid_io_fluid.WriteMesh(visualization_model_part.GetMesh());
            // gid_io_fluid.FinalizeMesh();
            // gid_io_fluid.InitializeResults(0, visualization_model_part.GetMesh());
            // gid_io_fluid.WriteNodalResults(DISTANCE, visualization_model_part.Nodes(), 0, 0);
            // gid_io_fluid.WriteNodalResults(VELOCITY, visualization_model_part.Nodes(), 0, 0);
            // gid_io_fluid.WriteNodalResults(PRESSURE, visualization_model_part.Nodes(), 0, 0);
            // gid_io_fluid.FinalizeResults();
        }

	    /** 
	     * Checks the embedded skin visualization process for tetrahedral meshes.
	     */
	    KRATOS_TEST_CASE_IN_SUITE(EmbeddedSkinVisualizationProcessTetrahedra, FluidDynamicsApplicationFastSuite)
		{

            // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
            Node<3>::Pointer p_point1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
            Node<3>::Pointer p_point2 = Kratos::make_shared<Node<3>>(2, 10.00, 0.00, 0.00);
            Node<3>::Pointer p_point3 = Kratos::make_shared<Node<3>>(3, 10.00, 10.00, 0.00);
            Node<3>::Pointer p_point4 = Kratos::make_shared<Node<3>>(4, 0.00, 10.00, 0.00);
            Node<3>::Pointer p_point5 = Kratos::make_shared<Node<3>>(5, 0.00, 0.00, 10.00);
            Node<3>::Pointer p_point6 = Kratos::make_shared<Node<3>>(6, 10.00, 0.00, 10.00);
            Node<3>::Pointer p_point7 = Kratos::make_shared<Node<3>>(7, 10.00, 10.00, 10.00);
            Node<3>::Pointer p_point8 = Kratos::make_shared<Node<3>>(8, 0.00, 10.00, 10.00);

            Hexahedra3D8<Node<3>> cube_domain(p_point1, p_point2, p_point3, p_point4, p_point5, p_point6, p_point7, p_point8);

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions" :  5,
                "element_name"        : "Element3D4N"
            })");

            ModelPart main_model_part("MainModelPart");
            main_model_part.AddNodalSolutionStepVariable(DISTANCE);
            main_model_part.AddNodalSolutionStepVariable(VELOCITY);
            main_model_part.AddNodalSolutionStepVariable(PRESSURE);
            StructuredMeshGeneratorProcess(cube_domain, main_model_part, mesher_parameters).Execute();

            // Set the nodal values
            const double level_set_height = 4.5;
            for (unsigned int i_node = 0; i_node < main_model_part.NumberOfNodes(); ++i_node){
                auto it_node = main_model_part.NodesBegin() + i_node;

                // Set DISTANCE values
                const double node_distance = (it_node->Z()) - level_set_height;
                it_node->FastGetSolutionStepValue(DISTANCE) = node_distance;

                // Set the VELOCITY and PRESSURE values
                double p_node;
                array_1d<double, 3> v_node = ZeroVector(3);
                if (node_distance > 0.0){
                    v_node[0] = it_node->X()*it_node->Z();
                    v_node[1] = it_node->X()*it_node->Z();
                    p_node = (std::pow(it_node->X(), 3))*it_node->Z();
                } else {
                    v_node[0] = -(it_node->X()*it_node->Z());
                    v_node[1] = -(it_node->X()*it_node->Z());
                    p_node = -(std::pow(it_node->X(), 2.5))*it_node->Z();
                }

                it_node->FastGetSolutionStepValue(VELOCITY) = v_node;
                it_node->FastGetSolutionStepValue(PRESSURE) = p_node;

            }

            // Create the visualization model part
            ModelPart visualization_model_part("VisualizationModelPart");
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

            // GidIO<> gid_io_fluid("/home/rzorrilla/Kratos/tests/visualizaton_model_part_3d", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
            // gid_io_fluid.InitializeMesh(0.0);
            // gid_io_fluid.WriteMesh(visualization_model_part.GetMesh());
            // gid_io_fluid.FinalizeMesh();
            // gid_io_fluid.InitializeResults(0, visualization_model_part.GetMesh());
            // gid_io_fluid.WriteNodalResults(DISTANCE, visualization_model_part.Nodes(), 0, 0);
            // gid_io_fluid.WriteNodalResults(VELOCITY, visualization_model_part.Nodes(), 0, 0);
            // gid_io_fluid.WriteNodalResults(PRESSURE, visualization_model_part.Nodes(), 0, 0);
            // gid_io_fluid.FinalizeResults();

	    }
    } // namespace Testing
}  // namespace Kratos.
