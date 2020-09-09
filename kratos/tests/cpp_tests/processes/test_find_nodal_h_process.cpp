//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/find_nodal_h_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

//         void GiDIODebugNodalH(ModelPart& ThisModelPart)
//         {
//             GidIO<> gid_io("TEST_GRADIENT", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(ThisModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, ThisModelPart.GetMesh());
//             gid_io.WriteNodalResults(NODAL_H, ThisModelPart.Nodes(), label, 0);
//         }

        /**
        * Checks the correct work of the nodal H compute
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalH1, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(this_model_part);

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsHistoricalVariable>(this_model_part);
            process.Execute();

//             // DEBUG
//             GiDIODebugNodalH(this_model_part);

            const double tolerance = 1.0e-4;
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(6)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
        }

        /**
        * Checks the correct work of the nodal H non-historical compute
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalH1NonHistorical, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(1);

            CppTestsUtilities::Create2DGeometry(this_model_part);

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(this_model_part);
            process.Execute();

            // // DEBUG
            // GiDIODebugNodalH(this_model_part);

            const double tolerance = 1.0e-4;
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(1)->GetValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(2)->GetValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(5)->GetValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(6)->GetValue(NODAL_H) - 1.0, tolerance);
        }

        /**
        * Checks the correct work of the nodal H compute
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalH2, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create3DGeometry(this_model_part);

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsHistoricalVariable>(this_model_part);
            process.Execute();

//             // DEBUG
//             GiDIODebugNodalH(this_model_part);

            const double tolerance = 1.0e-4;
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(5)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(9)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(10)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(11)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(this_model_part.pGetNode(12)->FastGetSolutionStepValue(NODAL_H) - 1.0, tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
