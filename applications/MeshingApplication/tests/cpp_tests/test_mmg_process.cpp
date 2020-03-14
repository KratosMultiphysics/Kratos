//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:   BSD License
//      Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
#ifdef INCLUDE_MMG
// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "containers/model.h"
#include "meshing_application_variables.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "custom_processes/mmg_process.h"
#include "processes/find_nodal_h_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

//         void GiDIODebugMMG(ModelPart& ThisModelPart)
//         {
//             GidIO<> gid_io("TEST_MMG", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(ThisModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, ThisModelPart.GetMesh());
//             gid_io.WriteNodalResultsNonHistorical(NODAL_H, ThisModelPart.Nodes(), label);
//             gid_io.WriteNodalFlags(ACTIVE, "ACTIVE", ThisModelPart.Nodes(), label);
//         }

        /**
        * Checks the correct work of the level set MMG process
        * Test triangle
        */

        KRATOS_TEST_CASE_IN_SUITE(MMGProcess1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

            // We set the flag to check that is transfered
            for (auto& i_elem : r_model_part.Elements())
                i_elem.Set(ACTIVE, true);

            // Set DISTANCE and other variables
            Vector ref_metric(3);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 0.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_2D, ref_metric);
            }

            // Compute remesh
            Parameters params = Parameters(R"({ "echo_level" : 0 })" );
            MmgProcess<MMGLibrary::MMG2D> mmg_process(r_model_part, params);
            mmg_process.Execute();

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
            process.Execute();

//             // DEBUG
//             GiDIODebugMMG(r_model_part);

            const double tolerance = 1.0e-4;
            for (auto& i_node : r_model_part.Nodes())
                if (i_node.X() < 0.001 || i_node.X() > 1.9999)
                    KRATOS_CHECK_LESS_EQUAL(i_node.GetValue(NODAL_H) - 1.0, tolerance);

            for (auto& i_elem : r_model_part.Elements())
                KRATOS_CHECK(i_elem.Is(ACTIVE));
        }

        /**
        * Checks the correct work of the level set MMG process
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(MMGProcess2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create3DGeometry(r_model_part, "Element3D4N");

            // We set the flag to check that is transfered
            for (auto& i_elem : r_model_part.Elements())
                i_elem.Set(ACTIVE, true);

            // Set DISTANCE and other variables
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_3D, ref_metric);
            }

            // Compute remesh
            Parameters params = Parameters(R"({ "echo_level" : 0 })" );
            MmgProcess<MMGLibrary::MMG3D> mmg_process(r_model_part, params);
            mmg_process.Execute();

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
            process.Execute();

//             // DEBUG
//             GiDIODebugMMG(r_model_part);

            double max = 0.0;
            for (auto& i_node : r_model_part.Nodes())
                if (i_node.GetValue(NODAL_H) > max)
                    max = i_node.GetValue(NODAL_H);

            const double tolerance = 1.0e-2;
            KRATOS_CHECK_LESS_EQUAL(std::abs(max - 1.0/std::sqrt(2.0))/max, tolerance);

            for (auto& i_elem : r_model_part.Elements())
                KRATOS_CHECK(i_elem.Is(ACTIVE));
        }
    } // namespace Testing
}  // namespace Kratos.
#endif
