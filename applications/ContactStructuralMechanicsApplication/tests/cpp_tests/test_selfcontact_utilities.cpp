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
#include "includes/gid_io.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/self_contact_utilities.h"

namespace Kratos
{
    namespace Testing
    {
        void GiDIOSelfContactDebug(ModelPart& rModelPart)
        {
            GidIO<> gid_io("TEST_SELFCONTACT_UTILITIES", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteConditionsOnly);
            const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(rModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, rModelPart.GetMesh());
            gid_io.WriteNodalFlags(MASTER, "MASTER", rModelPart.Nodes(), label);
            gid_io.WriteNodalFlags(SLAVE, "SLAVE", rModelPart.Nodes(), label);
            gid_io.WriteNodalResults(NORMAL, rModelPart.Nodes(), label, 0);
        }

        /**
         * This method can be used to create a 3D plane condition set
         */
        void SimpleCreateNewProblem3D(ModelPart& rModelPart)
        {
            // Creating nodes
            rModelPart.CreateNewNode(195, 4.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(196, 4.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(240, 4.0000E+00, 5.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(241, 5.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(242, 4.0000E+00, 4.7000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(243, 4.7000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(278, 5.0000E+00, 5.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(285, 4.7000E+00, 4.7000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(290, 4.0000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(292, 6.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(296, 6.0000E+00, 4.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(297, 4.0000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(332, 5.0000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(333, 6.0000E+00, 5.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(337, 4.7000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(338, 6.0000E+00, 4.7000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(371, 6.0000E+00, 6.0000E+00, 0.0000E+00);
            rModelPart.CreateNewNode(372, 6.0000E+00, 6.0000E+00, 0.0000E+00);

            // Creating properties
            Properties::Pointer p_cond_prop_0 = rModelPart.CreateNewProperties(0);
            Properties::Pointer p_cond_prop_1 = rModelPart.CreateNewProperties(1);

            // Creating conditions
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 1, {{297,242,285,337}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 2, {{337,285,338,372}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 3, {{242,196,243,285}}, p_cond_prop_0);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 4, {{285,243,296,338}}, p_cond_prop_0);

            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 5, {{371,333,278,332}}, p_cond_prop_1);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 6, {{332,278,240,290}}, p_cond_prop_1);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 7, {{333,292,241,278}}, p_cond_prop_1);
            rModelPart.CreateNewCondition("SurfaceCondition3D4N", 8, {{278,241,195,240}}, p_cond_prop_1);
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

            // Creating the pairs
            SelfContactUtilities::ComputeSelfContactPairing(r_model_part);
//             SelfContactUtilities::NotPredefinedMasterSlave(r_model_part);

            // DEBUG
            GiDIOSelfContactDebug(r_model_part);
        }

    } // namespace Testing
}  // namespace Kratos.
