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
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/checks.h"
#include "linear_solvers/linear_solver.h"
#include "processes/structured_mesh_generator_process.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "processes/calculate_embedded_nodal_variable_from_skin_process.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef CalculateEmbeddedNodalVariableFromSkinProcess<double, SparseSpaceType, LocalSpaceType, LinearSolverType> EmbeddedNodalVariableProcessDouble;
    typedef CalculateEmbeddedNodalVariableFromSkinProcess<array_1d<double, 3>, SparseSpaceType, LocalSpaceType, LinearSolverType> EmbeddedNodalVariableProcessArray;

    KRATOS_TEST_CASE_IN_SUITE(CalculateEmbeddedNodalVariableFromSkinProcessDouble, KratosCoreFastSuite)
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
        skin_part.GetNode(901).FastGetSolutionStepValue(TEMPERATURE) = 1.0;
        skin_part.GetNode(902).FastGetSolutionStepValue(TEMPERATURE) = 2.0;
        skin_part.GetNode(903).FastGetSolutionStepValue(TEMPERATURE) = 3.0;
        skin_part.GetNode(904).FastGetSolutionStepValue(TEMPERATURE) = 2.0;

        // Compute distance
        CalculateDistanceToSkinProcess<2>(surface_part, skin_part).Execute();

        // Compute the embedded nodal variable values
        EmbeddedNodalVariableProcessDouble emb_nod_var_from_skin_proc(
            surface_part,
            skin_part,
            TEMPERATURE,
            TEMPERATURE,
            "continuous");

        emb_nod_var_from_skin_proc.Execute();

        // Check values
        std::vector<std::size_t> check_nodes_ids = {19, 20, 28, 45, 46, 54};
        std::vector<double> expected_values = {0.883925, 1.42073, 1.69546, 2.57927, 3.11607, 3.47553};
        for (std::size_t i_node = 0; i_node < check_nodes_ids.size(); ++i_node) {
            const auto p_node = surface_part.pGetNode(check_nodes_ids[i_node]);
            KRATOS_CHECK_NEAR(p_node->FastGetSolutionStepValue(TEMPERATURE), expected_values[i_node], 1e-5);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(CalculateEmbeddedNodalVariableFromSkinProcessArray, KratosCoreFastSuite)
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
        surface_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        StructuredMeshGeneratorProcess(geometry, surface_part, mesher_parameters).Execute();

        // Generate the skin
        ModelPart &skin_part = current_model.CreateModelPart("Skin");
        skin_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        skin_part.CreateNewNode(901, 0.24, 0.34, 0.0);
        skin_part.CreateNewNode(902, 0.76, 0.34, 0.0);
        skin_part.CreateNewNode(903, 0.76, 0.66, 0.0);
        skin_part.CreateNewNode(904, 0.24, 0.66, 0.0);
        Properties::Pointer p_properties(new Properties(0));
        skin_part.CreateNewElement("Element2D2N", 901, {{901, 902}}, p_properties);
        skin_part.CreateNewElement("Element2D2N", 902, {{902, 903}}, p_properties);
        skin_part.CreateNewElement("Element2D2N", 903, {{903, 904}}, p_properties);
        skin_part.CreateNewElement("Element2D2N", 904, {{904, 901}}, p_properties);
        skin_part.GetNode(901).FastGetSolutionStepValue(DISPLACEMENT_X) = 1.0;
        skin_part.GetNode(902).FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.0;
        skin_part.GetNode(903).FastGetSolutionStepValue(DISPLACEMENT_X) = 3.0;
        skin_part.GetNode(904).FastGetSolutionStepValue(DISPLACEMENT_Y) = 2.0;

        // Compute distance
        CalculateDiscontinuousDistanceToSkinProcess<2>(surface_part, skin_part).Execute();

        // Compute the embedded nodal variable values
        EmbeddedNodalVariableProcessArray emb_nod_var_from_skin_proc(
            surface_part,
            skin_part,
            DISPLACEMENT,
            DISPLACEMENT,
            "discontinuous");

        emb_nod_var_from_skin_proc.Execute();

        // Check values
        const std::vector<std::size_t> check_nodes_ids = {19, 20, 28, 45, 46, 54};
        const std::vector<double> expected_values_x = {0.344868, 1.83755, 1.56283, 2.99609, 2.57702, 1.75274};
        const std::vector<double> expected_values_y = {0.539057, -0.416819, 0.132631, -0.416819, 0.539057, 1.72279};
        for (std::size_t i_node = 0; i_node < check_nodes_ids.size(); ++i_node) {
            const auto p_node = surface_part.pGetNode(check_nodes_ids[i_node]);
            KRATOS_CHECK_NEAR(p_node->FastGetSolutionStepValue(DISPLACEMENT_X), expected_values_x[i_node], 1e-5);
            KRATOS_CHECK_NEAR(p_node->FastGetSolutionStepValue(DISPLACEMENT_Y), expected_values_y[i_node], 1e-5);
        }
    }

}  // namespace Testing.
}  // namespace Kratos.
