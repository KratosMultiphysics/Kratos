//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
//
//

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "includes/gid_io.h"
#include "geometries/hexahedra_3d_8.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"

// TODO: REMOVE THAT ONCE THE IMPLEMENTATION IS FINISHED.
// TODO: BESIDES, MOVE THE TEST TO THE CORE 
#include "custom_processes/embedded_skin_visualization_process.h"

namespace Kratos {
namespace Testing {
    KRATOS_TEST_CASE_IN_SUITE(CubeInCubeDiscontinuousDistanceProcess, KratosCoreFastSuite)
    {
        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, -0.5, -0.5, -0.5);
        Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2,  0.5, -0.5, -0.5);
        Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3,  0.5,  0.5, -0.5);
        Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, -0.5,  0.5, -0.5);
        Node<3>::Pointer p_point_5 = Kratos::make_shared<Node<3>>(5, -0.5, -0.5,  0.5);
        Node<3>::Pointer p_point_6 = Kratos::make_shared<Node<3>>(6,  0.5, -0.5,  0.5);
        Node<3>::Pointer p_point_7 = Kratos::make_shared<Node<3>>(7,  0.5,  0.5,  0.5);
        Node<3>::Pointer p_point_8 = Kratos::make_shared<Node<3>>(8, -0.5,  0.5,  0.5);

        Hexahedra3D8<Node<3> > geometry(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions":   22,
            "element_name":     "Element3D4N"
        })");

        ModelPart volume_part("Volume");
        volume_part.AddNodalSolutionStepVariable(VELOCITY);
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

        // Generate the cube skin
        const double cube_radious = 0.25;
        ModelPart skin_part("Skin");
        skin_part.AddNodalSolutionStepVariable(VELOCITY);
        skin_part.CreateNewNode(1, -cube_radious, -cube_radious, -cube_radious);
        skin_part.CreateNewNode(2,  cube_radious, -cube_radious, -cube_radious);
        skin_part.CreateNewNode(3,  cube_radious,  cube_radious, -cube_radious);
        skin_part.CreateNewNode(4, -cube_radious,  cube_radious, -cube_radious);
        skin_part.CreateNewNode(5, -cube_radious, -cube_radious,  cube_radious);
        skin_part.CreateNewNode(6,  cube_radious, -cube_radious,  cube_radious);
        skin_part.CreateNewNode(7,  cube_radious,  cube_radious,  cube_radious);
        skin_part.CreateNewNode(8, -cube_radious,  cube_radious,  cube_radious);
        Properties::Pointer p_properties(new Properties(0));
        skin_part.CreateNewElement("Element3D3N",  1, { 1,2,3 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  2, { 1,3,4 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  3, { 5,6,7 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  4, { 5,7,8 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  5, { 3,6,2 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  6, { 3,7,6 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  7, { 4,5,1 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  8, { 4,8,5 }, p_properties);
        skin_part.CreateNewElement("Element3D3N",  9, { 3,4,8 }, p_properties);
        skin_part.CreateNewElement("Element3D3N", 10, { 3,8,7 }, p_properties);
        skin_part.CreateNewElement("Element3D3N", 11, { 2,1,5 }, p_properties);
        skin_part.CreateNewElement("Element3D3N", 12, { 2,5,6 }, p_properties);

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Save the analytic distance value for comparison
        for (unsigned int i = 0; i < volume_part.NumberOfNodes(); ++i) {
            auto it_node = volume_part.NodesBegin() + i;
            // Analytical cube distance
            const double d_yz = std::abs(it_node->X()) - cube_radious;
            const double d_xz = std::abs(it_node->Y()) - cube_radious;
            const double d_xy = std::abs(it_node->Z()) - cube_radious;
            double anal_dist = std::max(d_yz, d_xz);
            anal_dist = std::max(anal_dist, d_xy);
            it_node->FastGetSolutionStepValue(DISTANCE) = anal_dist;
        }

        // GidIO<> gid_io_fluid_analytic("/home/rzorrilla/Desktop/CubeInCubeDiscontinuousDistanceProcessAnalytic", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        // gid_io_fluid_analytic.InitializeMesh(0.00);
        // gid_io_fluid_analytic.WriteMesh(volume_part.GetMesh());
        // gid_io_fluid_analytic.FinalizeMesh();
        // gid_io_fluid_analytic.InitializeResults(0, volume_part.GetMesh());
        // gid_io_fluid_analytic.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
        // gid_io_fluid_analytic.FinalizeResults();

        // Save the minimum nodal value for visualization
        const unsigned int n_elems = volume_part.NumberOfElements();
        // for (unsigned int i = 0; i < n_elems; ++i) {
        //     auto it_elem = volume_part.ElementsBegin() + i;
        //     auto elemental_distance = it_elem->GetValue(ELEMENTAL_DISTANCES);
        //     auto &r_geom = it_elem->GetGeometry();

        //     for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j) {
        //         double &dist = r_geom[j].FastGetSolutionStepValue(DISTANCE);
        //         if (std::abs(dist) > std::abs(elemental_distance[j])) {
        //             if ((std::abs(r_geom[j].X()) < cube_radious) && (std::abs(r_geom[j].Y()) < cube_radious) && (std::abs(r_geom[j].Z()) < cube_radious)){
        //                 dist = -std::abs(elemental_distance[j]);
        //             } else {
        //                 dist = std::abs(elemental_distance[j]);                        
        //             }
        //         }
        //     }
        // }

        // GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/CubeInCubeDiscontinuousDistanceProcess", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        // gid_io_fluid.InitializeMesh(0.00);
        // gid_io_fluid.WriteMesh(volume_part.GetMesh());
        // gid_io_fluid.FinalizeMesh();
        // gid_io_fluid.InitializeResults(0, volume_part.GetMesh());
        // gid_io_fluid.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
        // gid_io_fluid.FinalizeResults();

        // ModelPart visualization_model_part("VolumeVisualization");
        // visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
        // visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
        // visualization_model_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        // Parameters visualization_parameters( R"(
        // {
        //     "shape_functions"                     : "ausas",
        //     "reform_model_part_at_each_time_step" : false,
        //     "visualization_variables"             : ["DISTANCE"]
        // })");
        // EmbeddedSkinVisualizationProcess visualization_tool(volume_part, visualization_model_part, visualization_parameters);
        // visualization_tool.ExecuteInitialize();
        // visualization_tool.ExecuteInitialize();
        // visualization_tool.ExecuteBeforeSolutionLoop();
        // visualization_tool.ExecuteInitializeSolutionStep();
        // visualization_tool.ExecuteBeforeOutputStep();
        // visualization_tool.ExecuteFinalizeSolutionStep();

        // GidIO<> gid_io_fluid_visualization("/home/rzorrilla/Desktop/CubeInCubeDiscontinuousDistanceProcessVisualization", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        // gid_io_fluid_visualization.InitializeMesh(0.00);
        // gid_io_fluid_visualization.WriteMesh(visualization_model_part.GetMesh());
        // gid_io_fluid_visualization.FinalizeMesh();
        // gid_io_fluid_visualization.InitializeResults(0, visualization_model_part.GetMesh());
        // gid_io_fluid_visualization.WriteNodalResults(DISTANCE, visualization_model_part.Nodes(), 0, 0);
        // gid_io_fluid_visualization.FinalizeResults();

        auto elems_begin = volume_part.ElementsBegin();
        for (unsigned int i = 0; i < n_elems; ++i) {
            auto it_elem = elems_begin + i;
            auto &r_geom = it_elem->GetGeometry();
            // Get the discontinuous distance values
            Vector elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
            // Check if the element is intersected
            bool is_intersected = false;
            for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j) {
                if (elem_dist[j] != 1.0) {
                    is_intersected = true;
                    break;
                }
            }
            // If the element is intersected, compute the analytical distance values and check
            if (is_intersected) {
                for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j) {
                    // Analytical cube distance
                    const double d_yz = std::abs(r_geom[j].X()) - cube_radious;
                    const double d_xz = std::abs(r_geom[j].Y()) - cube_radious;
                    const double d_xy = std::abs(r_geom[j].Z()) - cube_radious;
                    double anal_dist = std::max(d_yz, d_xz);
                    anal_dist = std::max(anal_dist, d_xy);
                    // Check analytical values against discontinuous distance ones
                    KRATOS_CHECK_NEAR(elem_dist[j], anal_dist, 1e-6);
                }
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(HorizontalPlaneDiscontinuousDistanceProcess, KratosCoreFastSuite)
    {
        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, -0.5, -0.5, -0.5);
        Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.5, -0.5, -0.5);
        Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 0.5, 0.5, -0.5);
        Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, -0.5, 0.5, -0.5);
        Node<3>::Pointer p_point_5 = Kratos::make_shared<Node<3>>(5, -0.5, -0.5, 0.5);
        Node<3>::Pointer p_point_6 = Kratos::make_shared<Node<3>>(6, 0.5, -0.5, 0.5);
        Node<3>::Pointer p_point_7 = Kratos::make_shared<Node<3>>(7, 0.5, 0.5, 0.5);
        Node<3>::Pointer p_point_8 = Kratos::make_shared<Node<3>>(8, -0.5, 0.5, 0.5);

        Hexahedra3D8<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions":   7,
            "element_name":     "Element3D4N"
        })");

        ModelPart volume_part("Volume");
        volume_part.AddNodalSolutionStepVariable(VELOCITY);
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        StructuredMeshGeneratorProcess(geometry, volume_part, mesher_parameters).Execute();

        // Generate the cube skin
        const double plane_height = 0.0;
        ModelPart skin_part("Skin");
        skin_part.AddNodalSolutionStepVariable(VELOCITY);
        skin_part.CreateNewNode(1, -1.0, -1.0, plane_height);
        skin_part.CreateNewNode(2,  1.0, -1.0, plane_height);
        skin_part.CreateNewNode(3,  1.0,  1.0, plane_height);
        skin_part.CreateNewNode(4, -1.0,  1.0, plane_height);
        Properties::Pointer p_properties(new Properties(0));
        skin_part.CreateNewElement("Element3D3N", 1, {1, 2, 4}, p_properties);
        skin_part.CreateNewElement("Element3D3N", 2, {2, 3, 4}, p_properties);

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Save the analytic distance value for comparison
        for (unsigned int i = 0; i < volume_part.NumberOfNodes(); ++i) {
            auto it_node = volume_part.NodesBegin() + i;
            it_node->FastGetSolutionStepValue(DISTANCE) = it_node->Z() - plane_height;
        }

        // GidIO<> gid_io_fluid_analytic("/home/rzorrilla/Desktop/HorizontalPlaneDiscontinuousDistanceProcessAnalytic", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        // gid_io_fluid_analytic.InitializeMesh(0.00);
        // gid_io_fluid_analytic.WriteMesh(volume_part.GetMesh());
        // gid_io_fluid_analytic.FinalizeMesh();
        // gid_io_fluid_analytic.InitializeResults(0, volume_part.GetMesh());
        // gid_io_fluid_analytic.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
        // gid_io_fluid_analytic.FinalizeResults();

        // auto n_elems = volume_part.NumberOfElements();
        // auto elems_begin = volume_part.ElementsBegin();
        // for (unsigned int i = 0; i < n_elems; ++i){
        //     auto it_elem = elems_begin + i;
        //     if (it_elem->Id() == 957)
        //         KRATOS_WATCH(*it_elem)
        //     auto &r_geom = it_elem->GetGeometry();
        //     auto &elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
        //     for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j){
        //         if(std::abs(elem_dist[j] - (r_geom[j].Z() - plane_height)) > 1e-6){
        //             if (it_elem->Is(TO_SPLIT))
        //                 std::cout << "Element: " << it_elem->Id() << " Node: " << j << " elem_dist[j]: " << elem_dist[j] << " anal_dist: " << r_geom[j].Z() - plane_height << std::endl;
        //         } 
        //     }
        // }

        // ModelPart visualization_model_part("VolumeVisualization");
        // visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
        // visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
        // visualization_model_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        // Parameters visualization_parameters(R"(
        // {
        //     "shape_functions"                     : "ausas",
        //     "reform_model_part_at_each_time_step" : false,
        //     "visualization_variables"             : ["DISTANCE"]
        // })");
        // EmbeddedSkinVisualizationProcess visualization_tool(volume_part, visualization_model_part, visualization_parameters);
        // visualization_tool.ExecuteInitialize();
        // visualization_tool.ExecuteInitialize();
        // visualization_tool.ExecuteBeforeSolutionLoop();
        // visualization_tool.ExecuteInitializeSolutionStep();
        // visualization_tool.ExecuteBeforeOutputStep();
        // visualization_tool.ExecuteFinalizeSolutionStep();

        // GidIO<> gid_io_fluid_visualization("/home/rzorrilla/Desktop/HorizontalPlaneDiscontinuousDistanceProcessVisualization", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        // gid_io_fluid_visualization.InitializeMesh(0.00);
        // gid_io_fluid_visualization.WriteMesh(visualization_model_part.GetMesh());
        // gid_io_fluid_visualization.FinalizeMesh();
        // gid_io_fluid_visualization.InitializeResults(0, visualization_model_part.GetMesh());
        // gid_io_fluid_visualization.WriteNodalResults(DISTANCE, visualization_model_part.Nodes(), 0, 0);
        // gid_io_fluid_visualization.FinalizeResults();

        // Results check
        auto n_elems = volume_part.NumberOfElements();
        auto elems_begin = volume_part.ElementsBegin();
        for (unsigned int i = 0; i < n_elems; ++i){
            auto it_elem = elems_begin + i;
            auto &r_geom = it_elem->GetGeometry();
            // Get the discontinuous distance values
            Vector elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
            // Check if the element is intersected
            bool is_intersected = false;
            for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j){
                if (elem_dist[j] != 1.0){
                    is_intersected = true;
                    break;
                }
            }
            // If the element is intersected, compute the analytical distance values and check
            if (is_intersected){
                for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j){
                    KRATOS_CHECK_NEAR(elem_dist[j], r_geom[j].Z() - plane_height, 1e-6);
                }
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(Element957DiscontinuousDistanceProcess, KratosCoreFastSuite)
    {
        // Generate the evil element
        ModelPart volume_part("Volume");
        volume_part.AddNodalSolutionStepVariable(VELOCITY);
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        volume_part.CreateNewNode(1, 0.214286, -0.357143, 0.0714286);
        volume_part.CreateNewNode(2, 0.214286, -0.214286, 0.0714286);
        volume_part.CreateNewNode(3, 0.357143, -0.214286, 0.0714286);
        volume_part.CreateNewNode(4, 0.214286, -0.357143, -0.0714286);
        Properties::Pointer p_properties_0(new Properties(0));
        volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

        // Generate the cube skin
        const double plane_height = 0.0;
        ModelPart skin_part("Skin");
        skin_part.AddNodalSolutionStepVariable(VELOCITY);
        skin_part.CreateNewNode(1, -1.0, -1.0, plane_height);
        skin_part.CreateNewNode(2,  1.0, -1.0, plane_height);
        skin_part.CreateNewNode(3,  1.0,  1.0, plane_height);
        skin_part.CreateNewNode(4, -1.0,  1.0, plane_height);
        Properties::Pointer p_properties_1(new Properties(1));
        skin_part.CreateNewElement("Element3D3N", 1, {1, 2, 4}, p_properties_1);
        skin_part.CreateNewElement("Element3D3N", 2, {2, 3, 4}, p_properties_1);

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Save the analytic distance value for comparison
        for (unsigned int i = 0; i < volume_part.NumberOfNodes(); ++i) {
            auto it_node = volume_part.NodesBegin() + i;
            it_node->FastGetSolutionStepValue(DISTANCE) = it_node->Z() - plane_height;
        }

        // GidIO<> gid_io_fluid_analytic("/home/rzorrilla/Desktop/Element957DiscontinuousDistanceProcessAnalytic", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        // gid_io_fluid_analytic.InitializeMesh(0.00);
        // gid_io_fluid_analytic.WriteMesh(volume_part.GetMesh());
        // gid_io_fluid_analytic.FinalizeMesh();
        // gid_io_fluid_analytic.InitializeResults(0, volume_part.GetMesh());
        // gid_io_fluid_analytic.WriteNodalResults(DISTANCE, volume_part.Nodes(), 0, 0);
        // gid_io_fluid_analytic.FinalizeResults();

        // ModelPart visualization_model_part("VolumeVisualization");
        // visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
        // visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
        // visualization_model_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        // Parameters visualization_parameters(R"(
        // {
        //     "shape_functions"                     : "ausas",
        //     "reform_model_part_at_each_time_step" : false,
        //     "visualization_variables"             : ["DISTANCE"]
        // })");
        // EmbeddedSkinVisualizationProcess visualization_tool(volume_part, visualization_model_part, visualization_parameters);
        // visualization_tool.ExecuteInitialize();
        // visualization_tool.ExecuteInitialize();
        // visualization_tool.ExecuteBeforeSolutionLoop();
        // visualization_tool.ExecuteInitializeSolutionStep();
        // visualization_tool.ExecuteBeforeOutputStep();
        // visualization_tool.ExecuteFinalizeSolutionStep();

        // GidIO<> gid_io_fluid_visualization("/home/rzorrilla/Desktop/Element957DiscontinuousDistanceProcessVisualization", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        // gid_io_fluid_visualization.InitializeMesh(0.00);
        // gid_io_fluid_visualization.WriteMesh(visualization_model_part.GetMesh());
        // gid_io_fluid_visualization.FinalizeMesh();
        // gid_io_fluid_visualization.InitializeResults(0, visualization_model_part.GetMesh());
        // gid_io_fluid_visualization.WriteNodalResults(DISTANCE, visualization_model_part.Nodes(), 0, 0);
        // gid_io_fluid_visualization.FinalizeResults();

        // Results check
        auto n_elems = volume_part.NumberOfElements();
        auto elems_begin = volume_part.ElementsBegin();
        for (unsigned int i = 0; i < n_elems; ++i){
            auto it_elem = elems_begin + i;
            auto &r_geom = it_elem->GetGeometry();
            // Get the discontinuous distance values
            Vector elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);
            // Check if the element is intersected
            bool is_intersected = false;
            for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j){
                if (elem_dist[j] != 1.0){
                    is_intersected = true;
                    break;
                }
            }
            // If the element is intersected, compute the analytical distance values and check
            if (is_intersected){
                for (unsigned int j = 0; j < r_geom.PointsNumber(); ++j){
                    KRATOS_CHECK_NEAR(elem_dist[j], r_geom[j].Z() - plane_height, 1e-6);
                }
            }
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximation2D, KratosCoreFastSuite)
    {
        // Generate the triangular element
        ModelPart volume_part("Volume");
        volume_part.AddNodalSolutionStepVariable(VELOCITY);
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        volume_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        volume_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        volume_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        Properties::Pointer p_properties_0(new Properties(0));
        volume_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);

        // Generate the skin such that one edge is cut twice
        ModelPart skin_part("Skin");
        skin_part.AddNodalSolutionStepVariable(VELOCITY);
        skin_part.CreateNewNode(1,  0.8, -0.1, 0.0);
        skin_part.CreateNewNode(2,  0.8,  0.8, 0.0);
        skin_part.CreateNewNode(3, -0.1,  0.8, 0.0);
        Properties::Pointer p_properties_1(new Properties(1));
        skin_part.CreateNewElement("Element2D2N", 1, {{1,2}}, p_properties_1);
        skin_part.CreateNewElement("Element2D2N", 2, {{2,3}}, p_properties_1);

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Call the visualization utility to see the resultant splitting pattern
        ModelPart visualization_model_part("VolumeVisualization");
        visualization_model_part.AddNodalSolutionStepVariable(VELOCITY);
        visualization_model_part.AddNodalSolutionStepVariable(DISTANCE);
        visualization_model_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        Parameters visualization_parameters(R"(
        {
            "shape_functions"                     : "ausas",
            "reform_model_part_at_each_time_step" : false,
            "visualization_variables"             : ["DISTANCE"]
        })");
        EmbeddedSkinVisualizationProcess visualization_tool(volume_part, visualization_model_part, visualization_parameters);
        visualization_tool.ExecuteInitialize();
        visualization_tool.ExecuteInitialize();
        visualization_tool.ExecuteBeforeSolutionLoop();
        visualization_tool.ExecuteInitializeSolutionStep();
        visualization_tool.ExecuteBeforeOutputStep();
        visualization_tool.ExecuteFinalizeSolutionStep();

        GidIO<> gid_io_fluid_visualization("/home/rzorrilla/Desktop/DiscontinuousDistanceProcessPlaneApproximation2D", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        gid_io_fluid_visualization.InitializeMesh(0.00);
        gid_io_fluid_visualization.WriteMesh(visualization_model_part.GetMesh());
        gid_io_fluid_visualization.FinalizeMesh();
        gid_io_fluid_visualization.InitializeResults(0, visualization_model_part.GetMesh());
        gid_io_fluid_visualization.WriteNodalResults(DISTANCE, visualization_model_part.Nodes(), 0, 0);
        gid_io_fluid_visualization.FinalizeResults();

    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximationSkewed, KratosCoreFastSuite)
    {
        // Generate the triangular element
        ModelPart volume_part("Volume");
        volume_part.AddNodalSolutionStepVariable(VELOCITY);
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        volume_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        volume_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        volume_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        volume_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        Properties::Pointer p_properties_0(new Properties(0));
        volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

        // Generate the skin such that one edge is cut twice
        ModelPart skin_part("Skin");
        skin_part.AddNodalSolutionStepVariable(VELOCITY);
        skin_part.CreateNewNode(1, -1.0, -1.0,  0.75);
        skin_part.CreateNewNode(2,  1.0, -1.0,  0.75);
        skin_part.CreateNewNode(3, -1.0,  1.0,  0.75);
        skin_part.CreateNewNode(4, 0.75, -1.0,  1.0);
        skin_part.CreateNewNode(5, 0.75, -1.0, -1.0);
        skin_part.CreateNewNode(6, 0.75,  1.0, -1.0);
        Properties::Pointer p_properties_1(new Properties(1));
        skin_part.CreateNewElement("Element3D3N", 1, {1,2,3}, p_properties_1);
        skin_part.CreateNewElement("Element3D3N", 2, {4,5,6}, p_properties_1);

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const Vector elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(elem_dist[0], -0.75, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[1], 0.25, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[2], -1.03078, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[3], 0.25, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessPlaneApproximationVertical, KratosCoreFastSuite)
    {
        // Generate the triangular element
        ModelPart volume_part("Volume");
        volume_part.AddNodalSolutionStepVariable(VELOCITY);
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY);
        volume_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        volume_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        volume_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        volume_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        Properties::Pointer p_properties_0(new Properties(0));
        volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

        // Generate the skin such that there is 4 intersection pts.
        // Recall that with more than 3 intersection pts. the plane 
        // approximation is used. Since the skin in here yields a 
        // uniplanar intersection, the approximated plane is the 
        // same one as the original intersection one.
        ModelPart skin_part("Skin");
        skin_part.AddNodalSolutionStepVariable(VELOCITY);
        skin_part.CreateNewNode(1, 0.5, -1.0,  1.0);
        skin_part.CreateNewNode(2, 0.5, -1.0, -1.0);
        skin_part.CreateNewNode(3, 0.5,  1.0, -1.0);
        skin_part.CreateNewNode(4, 0.5,  1.0,  1.0);
        Properties::Pointer p_properties_1(new Properties(1));
        skin_part.CreateNewElement("Element3D3N", 1, {1,2,4}, p_properties_1);
        skin_part.CreateNewElement("Element3D3N", 2, {2,3,4}, p_properties_1);

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Check values
        const Vector elem_dist = (volume_part.ElementsBegin())->GetValue(ELEMENTAL_DISTANCES);
        KRATOS_CHECK_NEAR(elem_dist[0], -0.5, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[1],  0.5, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[2], -0.5, 1e-10);
        KRATOS_CHECK_NEAR(elem_dist[3], -0.5, 1e-10);
    }

    KRATOS_TEST_CASE_IN_SUITE(DiscontinuousDistanceProcessOneEdgeIntersections, KratosCoreFastSuite)
    {
        // Generate the triangular element
        ModelPart volume_part("Volume");
        volume_part.AddNodalSolutionStepVariable(DISTANCE);
        volume_part.CreateNewNode(1, 0.666963, 0.800762, 0.388769);
        volume_part.CreateNewNode(2, 0.731067, 0.821936, 0.422077);
        volume_part.CreateNewNode(3, 0.652002, 0.85453, 0.463652);
        volume_part.CreateNewNode(4, 0.684484, 0.796908, 0.48275);
        Properties::Pointer p_properties_0(new Properties(0));
        volume_part.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_properties_0);

        // Generate the skin such that it only intersects in one edge
        ModelPart skin_part("Skin");
        skin_part.CreateNewNode(1, 0.675, 0.803109, 0.5);
        skin_part.CreateNewNode(2, 0.663088, 0.808771, 0.476277);
        skin_part.CreateNewNode(3, 0.685008, 0.796367, 0.479053);
        skin_part.CreateNewNode(4, 0.682845, 0.794215, 0.449949);
        Properties::Pointer p_properties_1(new Properties(1));
        skin_part.CreateNewElement("Element3D3N", 1, {1,2,3}, p_properties_1);
        skin_part.CreateNewElement("Element3D3N", 2, {3,2,4}, p_properties_1);

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess(volume_part, skin_part).Execute();

        // Call the visualization utility to see the resultant splitting pattern
        ModelPart visualization_model_part("VolumeVisualization");
        Parameters visualization_parameters(R"(
        {
            "shape_functions"                     : "ausas",
            "reform_model_part_at_each_time_step" : false,
            "visualization_variables"             : []
        })");
        EmbeddedSkinVisualizationProcess visualization_tool(volume_part, visualization_model_part, visualization_parameters);
        visualization_tool.ExecuteInitialize();
        visualization_tool.ExecuteInitialize();
        visualization_tool.ExecuteBeforeSolutionLoop();
        visualization_tool.ExecuteInitializeSolutionStep();
        visualization_tool.ExecuteBeforeOutputStep();
        visualization_tool.ExecuteFinalizeSolutionStep();

        GidIO<> gid_io_fluid_visualization("/home/rzorrilla/Desktop/DiscontinuousDistanceProcessOneEdgeIntersections", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        gid_io_fluid_visualization.InitializeMesh(0.00);
        gid_io_fluid_visualization.WriteMesh(visualization_model_part.GetMesh());
        gid_io_fluid_visualization.FinalizeMesh();
        gid_io_fluid_visualization.InitializeResults(0, visualization_model_part.GetMesh());
        gid_io_fluid_visualization.FinalizeResults();

    }

}  // namespace Testing.
}  // namespace Kratos.
