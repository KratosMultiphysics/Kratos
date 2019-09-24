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
#include "utilities/variable_utils.h"
#include "custom_utilities/self_contact_utilities.h"

namespace Kratos
{
    namespace Testing
    {
//         void GiDIOSelfContactDebug(ModelPart& rModelPart)
//         {
//             GidIO<> gid_io("TEST_SELFCONTACT_UTILITIES", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             gid_io.WriteNodalFlags(ACTIVE, "ACTIVE", rModelPart.Nodes(), label);
//             gid_io.WriteNodalFlags(MASTER, "MASTER", rModelPart.Nodes(), label);
//             gid_io.WriteNodalFlags(SLAVE, "SLAVE", rModelPart.Nodes(), label);
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
        * Checks the correct work of the weighted gap computation
        * Test 1
        */
        KRATOS_TEST_CASE_IN_SUITE(SelfContactUtilities1, KratosContactStructuralMechanicsFastSuite2)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Contact", 2);

            auto& r_process_info = r_model_part.GetProcessInfo();
            r_process_info[STEP] = 1;
            r_process_info[NL_ITERATION_NUMBER] = 1;

            // We create our problem
            SimpleCreateNewProblem3D(r_model_part);

            // All potential pairs
            for (auto& r_cond : r_model_part.Conditions()) {
                r_cond.SetValue(INDEX_MAP, Kratos::make_shared<IndexMap>());
            }
            for (auto& r_cond_1 : r_model_part.Conditions()) {
                auto p_pairs = r_cond_1.GetValue(INDEX_MAP);
                for (auto& r_cond_2 : r_model_part.Conditions()) {
                    if (r_cond_1.Id() != r_cond_2.Id()) {
                        p_pairs->AddId(r_cond_2.Id());
                    }
                }
            }

            // Creating the pairs
            SelfContactUtilities::ComputeSelfContactPairing(r_model_part);
//             SelfContactUtilities::NotPredefinedMasterSlave(r_model_part);

//             // DEBUG
//             GiDIOSelfContactDebug(r_model_part);

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

    } // namespace Testing
}  // namespace Kratos.
