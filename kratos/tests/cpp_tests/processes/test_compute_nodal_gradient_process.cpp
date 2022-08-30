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

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/cpp_tests_utilities.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"

/* Processes */
#include "processes/compute_nodal_gradient_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3>                                                    NodeType;

//         void GiDIODebugGradient(ModelPart& ThisModelPart)
//         {
//             GidIO<> gid_io("TEST_GRADIENT", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(ThisModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, ThisModelPart.GetMesh());
//             gid_io.WriteNodalResults(DISTANCE, ThisModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResults(DISTANCE_GRADIENT, ThisModelPart.Nodes(), label, 0);
//         }
        
        /**
        * Checks the correct work of the nodal gradient compute
        * Test triangle (2D)
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalGradient1, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(DISTANCE);
            this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create2DGeometry(this_model_part);

            // Initialize nodal area
            for (auto& node : this_model_part.Nodes())
                node.SetValue(NODAL_AREA, 0.0);

            // Set DISTANCE
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_AREA, 0.0);
            }

            typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> GradientType;
            GradientType process = GradientType(this_model_part, DISTANCE, DISTANCE_GRADIENT);
            process.Execute();

//             // DEBUG
//             GiDIODebugGradient(this_model_part);

            const double tolerance = 1.0e-8;
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(6)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
        }

        /**
        * Checks the correct work of the nodal gradient compute
        * Test triangle (3D)
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalGradient2, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(DISTANCE);
            this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create2DGeometry(this_model_part);

            // Initialize nodal area
            for (auto& node : this_model_part.Nodes())
                node.SetValue(NODAL_AREA, 0.0);

            // Set DISTANCE
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_AREA, 0.0);
            }

            typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> GradientType;
            GradientType process = GradientType(this_model_part, DISTANCE, DISTANCE_GRADIENT);
            process.Execute();

//             // DEBUG
//             GiDIODebugGradient(this_model_part);

            const double tolerance = 1.0e-8;
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(6)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
        }

        /**
        * Checks the correct work of the nodal gradient compute
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalGradient3, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(DISTANCE);
            this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create3DGeometry(this_model_part);

            // Initialize nodal area
            for (auto& node : this_model_part.Nodes())
                node.SetValue(NODAL_AREA, 0.0);

            // Set DISTANCE
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_AREA, 0.0);
            }

            // Compute gradient
            typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> GradientType;
            GradientType process = GradientType(this_model_part, DISTANCE, DISTANCE_GRADIENT);
            process.Execute();

//             // DEBUG
//             GiDIODebugGradient(this_model_part);

            const double tolerance = 1.0e-8;
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(3)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(9)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(10)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(11)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(12)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
        }

        /**
        * Checks the correct work of the nodal gradient compute
        * Test quadrilateral
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalGradient4, KratosCoreFastSuite)
        {
            Model current_model;
            
            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);
            
            this_model_part.AddNodalSolutionStepVariable(DISTANCE);
            this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);
            
            CppTestsUtilities::Create2DQuadrilateralsGeometry(this_model_part);
            
            // Initialize nodal area
            for (auto& node : this_model_part.Nodes())
                node.SetValue(NODAL_AREA, 0.0);
            
            // Set DISTANCE
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_AREA, 0.0);
            }
                         
            typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> GradientType;
            GradientType process = GradientType(this_model_part, DISTANCE, DISTANCE_GRADIENT);
            process.Execute();
            
//             // DEBUG         
//             GiDIODebugGradient(this_model_part);
            
            const double tolerance = 1.0e-8;
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(6)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
        }
        
        /** 
        * Checks the correct work of the nodal gradient compute
        * Test hexahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(NodalGradient5, KratosCoreFastSuite)
        {
            Model current_model;

            ModelPart& this_model_part = current_model.CreateModelPart("Main");
            this_model_part.SetBufferSize(2);

            this_model_part.AddNodalSolutionStepVariable(DISTANCE);
            this_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create3DHexahedraGeometry(this_model_part);

            // Initialize nodal area
            for (auto& node : this_model_part.Nodes())
                node.SetValue(NODAL_AREA, 0.0);

            // Set DISTANCE
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_AREA, 0.0);
            }

            // Compute gradient
            typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> GradientType;
            GradientType process = GradientType(this_model_part, DISTANCE, DISTANCE_GRADIENT);
            process.Execute();

//             // DEBUG
//             GiDIODebugGradient(this_model_part);

            const double tolerance = 1.0e-8;
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(1)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(2)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(3)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(5)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(9)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(10)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(11)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
            KRATOS_CHECK_LESS_EQUAL(std::abs(this_model_part.pGetNode(12)->FastGetSolutionStepValue(DISTANCE_GRADIENT_X)) - 1.0, tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
