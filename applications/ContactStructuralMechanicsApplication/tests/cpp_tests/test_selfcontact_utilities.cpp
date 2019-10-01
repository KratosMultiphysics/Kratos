// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:	   BSD License
//				   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "contact_structural_mechanics_application_variables.h"
// #include "includes/gid_io.h"
// #include "input_output/vtk_output.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/self_contact_utilities.h"

namespace Kratos
{
    namespace Testing
    {
//         void IOSelfContactDebug(ModelPart& rModelPart)
//         {
//             GidIO<> gid_io("TEST_SELFCONTACT_UTILITIES", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             gid_io.WriteNodalFlags(MASTER, "MASTER", rModelPart.Nodes(), label);
//             gid_io.WriteNodalFlags(SLAVE, "SLAVE", rModelPart.Nodes(), label);
//
//             Parameters io_parameters = Parameters(R"(
//             {
//                 "nodal_flags" : ["MASTER","SLAVE"]
//             })" );
//
//             VtkOutput(rModelPart, io_parameters).PrintOutput();
//         }

        /**
         * This method can be used to create a 3D plane condition set
         */
        void SimpleCreateNewProblem3D(ModelPart& rModelPart)
        {
            // Creating nodes
            rModelPart.CreateNewNode(1, 4.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(2, 4.0000E+00, 4.0000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(3, 4.0000E+00, 5.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(4, 5.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(5, 4.0000E+00, 4.7000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(6, 4.7000E+00, 4.0000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(7, 5.0000E+00, 5.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(8, 4.7000E+00, 4.7000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(9, 4.0000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(10, 6.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(11, 6.0000E+00, 4.0000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(12, 4.0000E+00, 6.0000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(13, 5.0000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(14, 6.0000E+00, 5.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(15, 4.7000E+00, 6.0000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(16, 6.0000E+00, 4.7000E+00, 0.5000E+00);
            rModelPart.CreateNewNode(17, 6.0000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(18, 6.0000E+00, 6.0000E+00, 0.5000E+00);

            // Creating properties
            Properties::Pointer p_cond_prop_0 = rModelPart.CreateNewProperties(0);
            Properties::Pointer p_cond_prop_1 = rModelPart.CreateNewProperties(1);

            // Creating conditions
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 1, {{12,5,8,15}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 2, {{15,8,16,18}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 3, {{5,2,6,8}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 4, {{8,6,11,16}}, p_cond_prop_0);

            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 5, {{17,14,7,13}}, p_cond_prop_1);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 6, {{13,7,3,9}}, p_cond_prop_1);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 7, {{14,10,4,7}}, p_cond_prop_1);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 8, {{7,4,1,3}}, p_cond_prop_1);
        }

        /**
         * This method can be used to create a 3D cylinder condition set
         */
        void CircularCreateNewProblem3D(ModelPart& rModelPart)
        {
            // Creating nodes
            rModelPart.CreateNewNode(1,-0.0097880904,-0.9999520955,0.5000000000);
            rModelPart.CreateNewNode(2,-0.0097880904,-0.9999520955,0.1000000000);
            rModelPart.CreateNewNode(3,-0.5223183884,-0.8527505503,0.5000000000);
            rModelPart.CreateNewNode(4,0.5082178134,-0.8612285725,0.5000000000);
            rModelPart.CreateNewNode(5,-0.5223183884,-0.8527505503,0.1000000000);
            rModelPart.CreateNewNode(6,0.5082178134,-0.8612285725,0.1000000000);
            rModelPart.CreateNewNode(7,-0.8832410687,-0.4689191983,0.5000000000);
            rModelPart.CreateNewNode(8,0.8770097528,-0.4804725731,0.5000000000);
            rModelPart.CreateNewNode(9,-0.8832410687,-0.4689191983,0.1000000000);
            rModelPart.CreateNewNode(10,0.8770097528,-0.4804725731,0.1000000000);
            rModelPart.CreateNewNode(11,-0.9992750647,0.0380702645,0.5000000000);
            rModelPart.CreateNewNode(12,0.9995947449,0.0284665770,0.5000000000);
            rModelPart.CreateNewNode(13,-0.9992750647,0.0380702645,0.1000000000);
            rModelPart.CreateNewNode(14,0.9995947449,0.0284665770,0.1000000000);
            rModelPart.CreateNewNode(15,-0.8489322851,0.5285016323,0.5000000000);
            rModelPart.CreateNewNode(16,0.8520420337,0.5234733735,0.5000000000);
            rModelPart.CreateNewNode(17,-0.8489322851,0.5285016323,0.1000000000);
            rModelPart.CreateNewNode(18,0.8520420337,0.5234733735,0.1000000000);
            rModelPart.CreateNewNode(19,-0.4820104280,0.8761654794,0.5000000000);
            rModelPart.CreateNewNode(20,0.4838042025,0.8751762643,0.5000000000);
            rModelPart.CreateNewNode(21,-0.4820104280,0.8761654794,0.1000000000);
            rModelPart.CreateNewNode(22,0.4838042025,0.8751762643,0.1000000000);
            rModelPart.CreateNewNode(23,0.0000000000,1.0000000000,0.5000000000);
            rModelPart.CreateNewNode(24,0.0000000000,1.0000000000,0.1000000000);

            // Creating properties
            Properties::Pointer p_cond_prop_0 = rModelPart.CreateNewProperties(0);

            // Creating conditions
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 1, {{23,24,21,19}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 2, {{19,21,17,15}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 3, {{15,17,13,11}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 4, {{11,13, 9, 7}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 5, {{7, 9, 5, 3}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 6, {{3, 5, 2, 1}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 7, {{1, 2, 6, 4}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 8, {{4, 6,10, 8}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 9, {{8,10,14,12}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 10, {{12,14,18,16}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 11, {{16,18,22,20}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 12, {{20,22,24,23}}, p_cond_prop_0);
        }

        /**
         * This method can be used to create a 3D an extended cylinder condition set
         */
        void ExtendedCircularCreateNewProblem3D(ModelPart& rModelPart)
        {
            // Creating nodes
            rModelPart.CreateNewNode( 1, 0.5877852508,-0.8090169955, 0.5000000000);
            rModelPart.CreateNewNode( 2, 0.5877852517,-0.8090169948, 0.2500000000);
            rModelPart.CreateNewNode( 3, 0.5877852508,-0.8090169955, 0.0000000000);
            rModelPart.CreateNewNode( 4,-0.0000000012,-1.0000000000, 0.5000000000);
            rModelPart.CreateNewNode( 5, 0.9510565162,-0.3090169946, 0.5000000000);
            rModelPart.CreateNewNode( 6, 0.0000000000,-1.0000000000, 0.2500000000);
            rModelPart.CreateNewNode( 7, 0.9510565166,-0.3090169936, 0.2500000000);
            rModelPart.CreateNewNode( 8,-0.0000000012,-1.0000000000, 0.0000000000);
            rModelPart.CreateNewNode( 9, 0.9510565162,-0.3090169946, 0.0000000000);
            rModelPart.CreateNewNode(10,-0.5877852527,-0.8090169941, 0.5000000000);
            rModelPart.CreateNewNode(11, 0.9510565166, 0.3090169935, 0.5000000000);
            rModelPart.CreateNewNode(12,-0.5877852517,-0.8090169948, 0.2500000000);
            rModelPart.CreateNewNode(13, 0.9510565162, 0.3090169946, 0.2500000000);
            rModelPart.CreateNewNode(14,-0.5877852527,-0.8090169941, 0.0000000000);
            rModelPart.CreateNewNode(15, 0.9510565166, 0.3090169935, 0.0000000000);
            rModelPart.CreateNewNode(16,-0.9510565169,-0.3090169926, 0.5000000000);
            rModelPart.CreateNewNode(17, 0.5877852536, 0.8090169934, 0.5000000000);
            rModelPart.CreateNewNode(18,-0.9510565166,-0.3090169936, 0.2500000000);
            rModelPart.CreateNewNode(19, 0.5877852527, 0.8090169941, 0.2500000000);
            rModelPart.CreateNewNode(20,-0.9510565169,-0.3090169926, 0.0000000000);
            rModelPart.CreateNewNode(21, 0.5877852536, 0.8090169934, 0.0000000000);
            rModelPart.CreateNewNode(22,-0.9510565159, 0.3090169957, 0.5000000000);
            rModelPart.CreateNewNode(23, 0.0000000000, 1.0000000000, 0.5000000000);
            rModelPart.CreateNewNode(24,-0.9510565162, 0.3090169946, 0.2500000000);
            rModelPart.CreateNewNode(25, 0.0000000000, 1.0000000000, 0.2500000000);
            rModelPart.CreateNewNode(26,-0.9510565159, 0.3090169957, 0.0000000000);
            rModelPart.CreateNewNode(27, 0.0000000000, 1.0000000000, 0.0000000000);
            rModelPart.CreateNewNode(28,-0.5877852517, 0.8090169948, 0.5000000000);
            rModelPart.CreateNewNode(29,-0.5877852527, 0.8090169941, 0.2500000000);
            rModelPart.CreateNewNode(30,-0.5877852517, 0.8090169948, 0.0000000000);

            // Creating properties
            Properties::Pointer p_cond_prop_0 = rModelPart.CreateNewProperties(0);

            // Creating conditions
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 1 , {{ 30 , 29 , 25 , 27 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 2 , {{ 26 , 24 , 29 , 30 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 3 , {{ 20 , 18 , 24 , 26 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 4 , {{ 14 , 12 , 18 , 20 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 5 , {{ 8 , 6 , 12 , 14 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 6 , {{ 3 , 2 , 6 , 8 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 7 , {{ 9 , 7 , 2 , 3 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 8 , {{ 15 , 13 , 7 , 9 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 9 , {{ 21 , 19 , 13 , 15 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 10 , {{ 27 , 25 , 19 , 21 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 11 , {{ 29 , 28 , 23 , 25 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 12 , {{ 24 , 22 , 28 , 29 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 13 , {{ 18 , 16 , 22 , 24 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 14 , {{ 12 , 10 , 16 , 18 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 15 , {{ 6 , 4 , 10 , 12 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 16 , {{ 2 , 1 , 4 , 6 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 17 , {{ 7 , 5 , 1 , 2 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 18 , {{ 13 , 11 , 5 , 7 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 19 , {{ 19 , 17 , 11 , 13 }}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 20 , {{ 25 , 23 , 17 , 19 }}, p_cond_prop_0);
        }

        /**
        * Checks the correct work of the self-contact utilities
        * Test 1
        */
        KRATOS_TEST_CASE_IN_SUITE(SelfContactUtilities1, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 2);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;

            // We create our problem
            SimpleCreateNewProblem3D(r_model_part);

            // All potential pairs
            SelfContactUtilities::FullAssignmentOfPairs(r_model_part);

            // Creating the pairs
            SelfContactUtilities::ComputeSelfContactPairing(r_model_part);

//             // DEBUG
//             IOSelfContactDebug(r_model_part);

            // Slave conditions
            KRATOS_CHECK(r_model_part.pGetCondition(1)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(1)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(3)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(3)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(4)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(4)->IsNot(MASTER));

            // Master conditions
            KRATOS_CHECK(r_model_part.pGetCondition(5)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(5)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->Is(MASTER));
        }

        /**
        * Checks the correct work of the self-contact utilities
        * Test 2
        */
        KRATOS_TEST_CASE_IN_SUITE(SelfContactUtilities2, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 2);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;

            // We create our problem
            CircularCreateNewProblem3D(r_model_part);

            // All potential pairs
            SelfContactUtilities::FullAssignmentOfPairs(r_model_part);

            // Creating the pairs
            SelfContactUtilities::ComputeSelfContactPairing(r_model_part);

//             // DEBUG
//             IOSelfContactDebug(r_model_part);

            // Slave conditions
            KRATOS_CHECK(r_model_part.pGetCondition(1)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(1)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(1)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(12)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(12)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(12)->IsNot(MASTER));

            // Inactive conditions
            KRATOS_CHECK(r_model_part.pGetCondition(3)->IsNot(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(3)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(3)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(11)->IsNot(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(11)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(11)->Is(MASTER));

            // Master conditions
            KRATOS_CHECK(r_model_part.pGetCondition(4)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(4)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(4)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(5)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(5)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(5)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(9)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(9)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(9)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(10)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(10)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(10)->Is(MASTER));
        }

        /**
        * Checks the correct work of the self-contact utilities
        * Test 3
        */
        KRATOS_TEST_CASE_IN_SUITE(SelfContactUtilities3, KratosContactStructuralMechanicsFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 2);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;

            // We create our problem
            ExtendedCircularCreateNewProblem3D(r_model_part);

            // All potential pairs
            SelfContactUtilities::FullAssignmentOfPairs(r_model_part);

            // Creating the pairs
            SelfContactUtilities::ComputeSelfContactPairing(r_model_part);

//             // DEBUG
//             IOSelfContactDebug(r_model_part);

            // Slave conditions
            KRATOS_CHECK(r_model_part.pGetCondition(1)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(1)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(1)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(11)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(11)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(11)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(2)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(12)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(12)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(12)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(10)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(10)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(10)->IsNot(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(20)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(20)->Is(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(20)->IsNot(MASTER));

            // Inactive conditions
            KRATOS_CHECK(r_model_part.pGetCondition(3)->IsNot(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(3)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(3)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(13)->IsNot(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(13)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(13)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(9)->IsNot(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(9)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(9)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(19)->IsNot(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(19)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(19)->Is(MASTER));

            // Master conditions
            KRATOS_CHECK(r_model_part.pGetCondition(4)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(4)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(4)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(14)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(14)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(14)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(5)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(5)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(5)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(15)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(15)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(15)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(6)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(16)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(16)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(16)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(7)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(17)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(17)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(17)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(8)->Is(MASTER));
            KRATOS_CHECK(r_model_part.pGetCondition(18)->Is(ACTIVE));
            KRATOS_CHECK(r_model_part.pGetCondition(18)->IsNot(SLAVE));
            KRATOS_CHECK(r_model_part.pGetCondition(18)->Is(MASTER));
        }

    } // namespace Testing
}  // namespace Kratos.
