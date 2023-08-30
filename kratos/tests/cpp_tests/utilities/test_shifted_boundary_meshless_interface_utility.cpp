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

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/generic_find_elements_neighbours_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/shifted_boundary_meshless_interface_utility.h"

namespace Kratos::Testing
{

namespace
{
    void AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest(
        Model& rModel,
        Parameters TestSettings)
    {
        // Generate a surface mesh (done with the StructuredMeshGeneratorProcess)
		auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
		auto p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
		auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
		auto p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
		Quadrilateral2D4<Node> square_geometry(p_point_1, p_point_2, p_point_3, p_point_4);

		Parameters mesher_parameters(R"(
		{
			"number_of_divisions":   7,
			"element_name": "Element2D3N",
            "create_skin_sub_model_part": false
		})");
		auto& r_test_model_part = rModel.CreateModelPart("TestModelPart");
        r_test_model_part.GetProcessInfo()[DOMAIN_SIZE] = 2;
        r_test_model_part.AddNodalSolutionStepVariable(DISTANCE);
		StructuredMeshGeneratorProcess(square_geometry, r_test_model_part, mesher_parameters).Execute();

        // Set distance field
        const double y_0 = 0.5;
        for (auto& r_node : r_test_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(DISTANCE) = r_node.Y() - y_0;
        }

        // Calculate the nodal neighbours
        // Note that this is required by the SBM WTE extension operator utility
        auto gl_nodal_neigh_proc = FindGlobalNodalNeighboursProcess(r_test_model_part);
        gl_nodal_neigh_proc.Execute();
        auto gen_elem_neigh_proc = GenericFindElementalNeighboursProcess(r_test_model_part);
        gen_elem_neigh_proc.Execute();

        // Call the utility calculating the SBM WTE extension operator
        auto p_sbm_wte_ext_op_utility = Kratos::make_unique<ShiftedBoundaryMeshlessInterfaceUtility>(rModel, TestSettings);
        p_sbm_wte_ext_op_utility->CalculateExtensionOperator();
    }
}

    KRATOS_TEST_CASE_IN_SUITE(ShiftedBoundaryMeshlessInterfaceUtilityConformingMLS, KratosCoreFastSuite)
    {
        // Test settings
        Model test_model;
        Parameters test_settings(R"({
            "model_part_name" : "TestModelPart",
            "boundary_sub_model_part_name" : "BoundaryModelPart",
            "sbm_interface_condition_name" : "GenericCondition",
            "conforming_basis" : true,
            "extension_operator_type" : "MLS"
        })");

        // Call the auxiliary test execution function
        AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest(test_model, test_settings);

        // Check results from the first SBM WTE condition
        const auto cond_begin = test_model.GetModelPart("TestModelPart.BoundaryModelPart").ConditionsBegin();
        const auto& r_N = cond_begin->GetValue(SHAPE_FUNCTIONS_VECTOR);
        const auto& r_DN_DX = cond_begin->GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);

        const double tolerance = 1.0e-8;
        const std::vector<double> expected_N({0.690217814756,0.103088070081,-0.033541375515,0.3223921909,0.0440525346534,-0.0722263529105,0.129993178184,-0.00796714246957,-0.0954187416683,0.0579677282616,-0.0403152864677,-0.0982426178045});
        Matrix expected_DN_DX(12,2);
        expected_DN_DX(0,0) = -7; expected_DN_DX(0,1) = 2.85767653557;
        expected_DN_DX(1,0) = 0; expected_DN_DX(1,1) = -1.44323298113;
        expected_DN_DX(2,0) = 0; expected_DN_DX(2,1) = 0.469579257211;
        expected_DN_DX(3,0) = 7; expected_DN_DX(3,1) = -3.03421661476;
        expected_DN_DX(4,0) = 0; expected_DN_DX(4,1) = -0.616735485148;
        expected_DN_DX(5,0) = 0; expected_DN_DX(5,1) = 1.01116894075;
        expected_DN_DX(6,0) = 0; expected_DN_DX(6,1) = -1.81990449457;
        expected_DN_DX(7,0) = 0; expected_DN_DX(7,1) = 0.111539994574;
        expected_DN_DX(8,0) = 0; expected_DN_DX(8,1) = 1.33586238336;
        expected_DN_DX(9,0) = 0; expected_DN_DX(9,1) = -0.811548195663;
        expected_DN_DX(10,0) = 0; expected_DN_DX(10,1) = 0.564414010548;
        expected_DN_DX(11,0) = 0; expected_DN_DX(11,1) = 1.37539664926;
        KRATOS_EXPECT_VECTOR_NEAR(r_N, expected_N, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_DN_DX, expected_DN_DX, tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(ShiftedBoundaryMeshlessInterfaceUtilityNonConformingMLS, KratosCoreFastSuite)
    {
        // Test settings
        Model test_model;
        Parameters test_settings(R"({
            "model_part_name" : "TestModelPart",
            "boundary_sub_model_part_name" : "BoundaryModelPart",
            "sbm_interface_condition_name" : "GenericCondition",
            "conforming_basis" : false,
            "extension_operator_type" : "MLS"
        })");

        // Call the auxiliary test execution function
        AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest(test_model, test_settings);

        // Check results from the first SBM WTE condition
        const auto cond_begin = test_model.GetModelPart("TestModelPart.BoundaryModelPart").ConditionsBegin();
        const auto& r_N = cond_begin->GetValue(SHAPE_FUNCTIONS_VECTOR);
        const auto& r_DN_DX = cond_begin->GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);

        const double tolerance = 1.0e-8;
        const std::vector<double> expected_N({0.850500942829,-0.0440476158159,0.495589016468,-0.214158103197,0.153910040703,-0.241794280987});
        Matrix expected_DN_DX(6,2);
        expected_DN_DX(0,0) = -4.16531166617; expected_DN_DX(0,1) = -3.11225503535;
        expected_DN_DX(1,0) = -0.422817410056; expected_DN_DX(1,1) = 3.20331597666;
        expected_DN_DX(2,0) = 1.8246012174; expected_DN_DX(2,1) = -2.4872855784;
        expected_DN_DX(3,0) = 0.351656935045; expected_DN_DX(3,1) = 2.30516369579;
        expected_DN_DX(4,0) = 2.34071044876; expected_DN_DX(4,1) = -1.40045938624;
        expected_DN_DX(5,0) = 0.071160475011; expected_DN_DX(5,1) = 1.49152032755;
        KRATOS_EXPECT_VECTOR_NEAR(r_N, expected_N, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_DN_DX, expected_DN_DX, tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(ShiftedBoundaryMeshlessInterfaceUtilityRBF, KratosCoreFastSuite)
    {
        // Test settings
        Model test_model;
        Parameters test_settings(R"({
            "model_part_name" : "TestModelPart",
            "boundary_sub_model_part_name" : "BoundaryModelPart",
            "sbm_interface_condition_name" : "GenericCondition",
            "conforming_basis" : true,
            "extension_operator_type" : "RBF"
        })");

        // Call the auxiliary test execution function
        AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest(test_model, test_settings);

        // Check results from the first SBM WTE condition
        const auto cond_begin = test_model.GetModelPart("TestModelPart.BoundaryModelPart").ConditionsBegin();
        const auto& r_N = cond_begin->GetValue(SHAPE_FUNCTIONS_VECTOR);
        const auto& r_DN_DX = cond_begin->GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);

        const double tolerance = 1.0e-8;
        const std::vector<double> expected_N({2.19812295036,-1.84616655368,0.555840635338,-0.747583442815,1.02822480988,-0.207155100513,0.600496011945,-0.725928425621,0.148970755181,-0.107282272913,0.129526737242,-0.0270661044029});
        Matrix expected_DN_DX(12,2);
        expected_DN_DX(0,0) = -7; expected_DN_DX(0,1) = -18.2529953628;
        expected_DN_DX(1,0) = 0; expected_DN_DX(1,1) = 25.8463317515;
        expected_DN_DX(2,0) = 0; expected_DN_DX(2,1) = -7.78176889473;
        expected_DN_DX(3,0) = 7; expected_DN_DX(3,1) = 11.9454422572;
        expected_DN_DX(4,0) = 0; expected_DN_DX(4,1) = -14.3951473383;
        expected_DN_DX(5,0) = 0; expected_DN_DX(5,1) = 2.90017140718;
        expected_DN_DX(6,0) = 0; expected_DN_DX(6,1) = -8.40694416722;
        expected_DN_DX(7,0) = 0; expected_DN_DX(7,1) = 10.1629979587;
        expected_DN_DX(8,0) = 0; expected_DN_DX(8,1) = -2.08559057253;
        expected_DN_DX(9,0) = 0; expected_DN_DX(9,1) = 1.50195182078;
        expected_DN_DX(10,0) = 0; expected_DN_DX(10,1) = -1.81337432138;
        expected_DN_DX(11,0) = 0; expected_DN_DX(11,1) = 0.378925461641;
        KRATOS_EXPECT_VECTOR_NEAR(r_N, expected_N, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_DN_DX, expected_DN_DX, tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(ShiftedBoundaryMeshlessInterfaceUtilityGradientBased, KratosCoreFastSuite)
    {
        // Test settings
        Model test_model;
        Parameters test_settings(R"({
            "model_part_name" : "TestModelPart",
            "boundary_sub_model_part_name" : "BoundaryModelPart",
            "sbm_interface_condition_name" : "GenericCondition",
            "conforming_basis" : true,
            "extension_operator_type" : "gradient_based"
        })");

        // Call the auxiliary test execution function
        AuxiliaryShiftedBoundaryMeshlessInterfaceUtilityTest(test_model, test_settings);

        // Check results from the first SBM WTE condition
        const auto cond_begin = test_model.GetModelPart("TestModelPart.BoundaryModelPart").ConditionsBegin();
        const auto& r_N = cond_begin->GetValue(SHAPE_FUNCTIONS_VECTOR);
        const auto& r_DN_DX = cond_begin->GetValue(SHAPE_FUNCTIONS_GRADIENT_MATRIX);

        const double tolerance = 1.0e-8;
        const std::vector<double> expected_N({0.902712989246,-0.146446609407,0.597287010754,-0.215482203136,2.29934717029e-17,-0.138071187458});
        Matrix expected_DN_DX(6,2);
        expected_DN_DX(0,0) = -7; expected_DN_DX(0,1) = -0.117255907286;
        expected_DN_DX(1,0) = 0; expected_DN_DX(1,1) = 2.05025253169;
        expected_DN_DX(2,0) = 7; expected_DN_DX(2,1) = -6.88274409271;
        expected_DN_DX(3,0) = 0; expected_DN_DX(3,1) = 3.0167508439;
        expected_DN_DX(4,0) = 0; expected_DN_DX(4,1) = -3.21908603841e-16;
        expected_DN_DX(5,0) = 0; expected_DN_DX(5,1) = 1.93299662441;
        KRATOS_EXPECT_VECTOR_NEAR(r_N, expected_N, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(r_DN_DX, expected_DN_DX, tolerance);
    }

} // namespace Kratos::Testing