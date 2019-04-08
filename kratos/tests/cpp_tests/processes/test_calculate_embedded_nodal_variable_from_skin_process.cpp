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

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/gid_io.h"
#include "linear_solvers/linear_solver.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "processes/calculate_embedded_nodal_variable_from_skin_process.h"
#include "spaces/ublas_space.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(CalculateEmbeddedNodalVariableFromSkinProcessHorizontalPlane2D, KratosCoreFastSuite)
    {
        Model current_model;

        // Generate the element
        ModelPart &fluid_part = current_model.CreateModelPart("Surface");
        fluid_part.AddNodalSolutionStepVariable(DISTANCE);
        fluid_part.AddNodalSolutionStepVariable(TEMPERATURE);
        fluid_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        fluid_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        fluid_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        fluid_part.CreateNewNode(4, 1.0, 1.0, 0.0);
        Properties::Pointer p_properties_0(new Properties(0));
        fluid_part.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_properties_0);
        fluid_part.CreateNewElement("Element2D3N", 2, {2, 4, 3}, p_properties_0);

        // Generate the skin
        const double plane_height = 0.5;
        ModelPart &skin_part = current_model.CreateModelPart("Skin");
        skin_part.AddNodalSolutionStepVariable(TEMPERATURE);
        skin_part.CreateNewNode(1, -1.0, plane_height, 0.0);
        skin_part.CreateNewNode(2,  0.6, plane_height, 0.0);
        Properties::Pointer p_properties_1(new Properties(1));
        skin_part.CreateNewElement("Element2D2N", 1, {{1, 2}}, p_properties_1);
        skin_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE) = 1.0;
        skin_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE) = 2.0;

        // Compute the discontinuous distance function
        CalculateDiscontinuousDistanceToSkinProcess<2> disc_dist_proc(fluid_part, skin_part);
        disc_dist_proc.Execute();

        // Compute the embedded nodal variable values
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef CalculateEmbeddedNodalVariableFromSkinProcess<2, double, SparseSpaceType, LocalSpaceType, LinearSolverType> EmbeddedNodalVariableProcess;
        EmbeddedNodalVariableProcess emb_nod_var_from_skin_proc(
            fluid_part,
            skin_part,
            TEMPERATURE,
            TEMPERATURE,
            "discontinuous");

        emb_nod_var_from_skin_proc.Execute();
        emb_nod_var_from_skin_proc.Clear();

        // Check values
        KRATOS_WATCH(fluid_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(fluid_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(fluid_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(fluid_part.GetNode(4).FastGetSolutionStepValue(TEMPERATURE))
        std::vector<double> expected_values = {0.0, 0.0, 0.0, 0.0};
        for (std::size_t i_node = 0; i_node < fluid_part.NumberOfNodes(); ++i_node) {
            const auto it_node = fluid_part.NodesBegin() + i_node;
            KRATOS_CHECK_NEAR(it_node->FastGetSolutionStepValue(TEMPERATURE), expected_values[i_node], 1e-12);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(CalculateEmbeddedNodalVariableFromSkinProcess2D, KratosCoreFastSuite)
    {
        // Generate a volume mesh (done with the StructuredMeshGeneratorProcess)
        Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.00, 0.00, 0.00);
        Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.00, 1.00, 0.00);
        Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 1.00, 1.00, 0.00);
        Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 1.00, 0.00, 0.00);

        Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

        Parameters mesher_parameters(R"({
			"number_of_divisions": 7,
			"element_name": "Element2D3N"
		})");

        Model current_model;
        ModelPart &surface_part = current_model.CreateModelPart("Volume");
        surface_part.AddNodalSolutionStepVariable(DISTANCE);
        surface_part.AddNodalSolutionStepVariable(TEMPERATURE);
        StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();

        // Generate the skin
        ModelPart &skin_part = current_model.CreateModelPart("Skin");
        skin_part.AddNodalSolutionStepVariable(TEMPERATURE);
        skin_part.CreateNewNode(901, 0.24, 0.34, 0.0);
        skin_part.CreateNewNode(902, 0.76, 0.34, 0.0);
        skin_part.CreateNewNode(903, 0.76, 0.66, 0.0);
        skin_part.CreateNewNode(904, 0.24, 0.66, 0.0);
        Properties::Pointer p_properties(new Properties(0));
        skin_part.CreateNewElement("Element2D2N", 901, {{901, 902}}, p_properties);
        skin_part.CreateNewElement("Element2D2N", 902, {{902, 903}}, p_properties);
        skin_part.CreateNewElement("Element2D2N", 903, {{903, 904}}, p_properties);
        skin_part.CreateNewElement("Element2D2N", 904, {{904, 901}}, p_properties);
        skin_part.GetNode(901).FastGetSolutionStepValue(TEMPERATURE) = 4.0;
        skin_part.GetNode(902).FastGetSolutionStepValue(TEMPERATURE) = 4.0;
        skin_part.GetNode(903).FastGetSolutionStepValue(TEMPERATURE) = 4.0;
        skin_part.GetNode(904).FastGetSolutionStepValue(TEMPERATURE) = 4.0;

        // Compute distance
        CalculateDistanceToSkinProcess<2>(surface_part, skin_part).Execute();

        // Compute the embedded nodal variable values
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef CalculateEmbeddedNodalVariableFromSkinProcess<2, double, SparseSpaceType, LocalSpaceType, LinearSolverType> EmbeddedNodalVariableProcess;
        EmbeddedNodalVariableProcess emb_nod_var_from_skin_proc(
            surface_part,
            skin_part,
            TEMPERATURE,
            TEMPERATURE,
            "continuous");

        emb_nod_var_from_skin_proc.Execute();
        emb_nod_var_from_skin_proc.Clear();

        GidIO<> gid_io_fluid("/home/rzorrilla/Desktop/surface_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        gid_io_fluid.InitializeMesh(0.00);
        gid_io_fluid.WriteMesh(surface_part.GetMesh());
        gid_io_fluid.FinalizeMesh();
        gid_io_fluid.InitializeResults(0, surface_part.GetMesh());
        gid_io_fluid.WriteNodalResults(DISTANCE, surface_part.Nodes(), 0, 0);
        gid_io_fluid.WriteNodalResults(TEMPERATURE, surface_part.Nodes(), 0, 0);
        gid_io_fluid.FinalizeResults();

        GidIO<> gid_io_skin("/home/rzorrilla/Desktop/skin_mesh", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
        gid_io_skin.InitializeMesh(0.00);
        gid_io_skin.WriteMesh(skin_part.GetMesh());
        gid_io_skin.FinalizeMesh();
        gid_io_skin.InitializeResults(0, skin_part.GetMesh());
        gid_io_skin.WriteNodalResults(TEMPERATURE, skin_part.Nodes(), 0, 0);
        gid_io_skin.FinalizeResults();

        // Check values
        KRATOS_WATCH(surface_part.GetNode(1).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(surface_part.GetNode(2).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(surface_part.GetNode(3).FastGetSolutionStepValue(TEMPERATURE))
        KRATOS_WATCH(surface_part.GetNode(4).FastGetSolutionStepValue(TEMPERATURE))
        std::vector<double> expected_values = {0.0, 0.0, 0.0, 0.0};
        for (std::size_t i_node = 0; i_node < surface_part.NumberOfNodes(); ++i_node) {
            const auto it_node = surface_part.NodesBegin() + i_node;
            KRATOS_CHECK_NEAR(it_node->FastGetSolutionStepValue(TEMPERATURE), expected_values[i_node], 1e-12);
        }
    }

}  // namespace Testing.
}  // namespace Kratos.
