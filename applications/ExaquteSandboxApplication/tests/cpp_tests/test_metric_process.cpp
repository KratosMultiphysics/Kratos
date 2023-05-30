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
//                   Riccardo Tosi
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "containers/model.h"

/* Processes */
#include "processes/compute_nodal_gradient_process.h"
#include "custom_processes/metrics_divergencefree_process.h"
#include "meshing_application_variables.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node NodeType;

//         void GiDIODebugMetric(ModelPart& rModelPart)
//         {
//             GidIO<> gid_io("TEST_METRIC", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             gid_io.WriteNodalResults(DISTANCE, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResults(DISTANCE_GRADIENT, rModelPart.Nodes(), label, 0);
//             gid_io.WriteNodalResultsNonHistorical(METRIC_TENSOR_2D, rModelPart.Nodes(), label);
// //             gid_io.WriteNodalResultsNonHistorical(METRIC_TENSOR_3D, rModelPart.Nodes(), label); // NOTE: 6 components not suported, update
//         }
//
//         void GiDIODebugMetricSPR(ModelPart& rModelPart)
//         {
//             GidIO<> gid_io("TEST_METRIC_SPR", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             gid_io.WriteNodalResults(DISPLACEMENT, rModelPart.Nodes(), label, 0);
//             gid_io.PrintOnGaussPoints(ERROR_INTEGRATION_POINT, rModelPart, label);
//             gid_io.PrintOnGaussPoints(CAUCHY_STRESS_VECTOR, rModelPart, label);
//             gid_io.PrintOnGaussPoints(STRAIN_ENERGY, rModelPart, label);
//             gid_io.WriteNodalResultsNonHistorical(METRIC_TENSOR_2D, rModelPart.Nodes(), label);
// //             gid_io.WriteNodalResultsNonHistorical(METRIC_TENSOR_3D, rModelPart.Nodes(), label); // NOTE: 6 components not suported, update
//         }

        void Create2DGeometry(ModelPart& rModelPart, const std::string& ElementName)
        {
            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            // First we create the nodes
            rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Now we create the elements
            rModelPart.CreateNewElement(ElementName, 1, {{1,2,3}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 2, {{1,3,4}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 3, {{2,5,3}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 4, {{5,6,3}}, p_elem_prop);
        }

        void Create3DGeometry(ModelPart& rModelPart, const std::string& ElementName)
        {
            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            // First we create the nodes
            rModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            rModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            rModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            rModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

            rModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            rModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            rModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            rModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

            // Now we create the elements
            rModelPart.CreateNewElement(ElementName, 1, {{12,10,8,9}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 2, {{4,6,9,7}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 3, {{11,7,9,8}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 4, {{5,3,8,6}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 5, {{4,6,7,3}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 6, {{2,3,5,6}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 7, {{10,9,6,8}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 8, {{7,8,3,6}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 9, {{7,8,6,9}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 10, {{4,1,6,3}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 11, {{9,12,11,8}}, p_elem_prop);
            rModelPart.CreateNewElement(ElementName, 12, {{3,2,1,6}}, p_elem_prop);
        }


        /**
        * Checks the correct work of the divergencefree metric process with maximum strategy
        * Test triangle
        */

        KRATOS_TEST_CASE_IN_SUITE(DivergenceFreeMetricProcess1, ExaquteSandboxApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            Parameters this_parameters = Parameters(R"(
            {
                "minimal_size"                        : 0.01,
                "maximal_size"                        : 10.0,
                "refinement_strategy"                 : "maximum_strategy",
                "reference_variable_name"             : "DISTANCE"
            })"
            );

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            Create2DGeometry(r_model_part, "Element2D3N");

            // Set DISTANCE variable
            for (std::size_t i_elem = 0; i_elem < r_model_part.Elements().size(); ++i_elem) {
                auto it_elem = r_model_part.Elements().begin() + i_elem;
                it_elem->SetValue(DISTANCE, 1.0);
            }

            // Set other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_2D, ZeroVector(3));
            }

            // Compute metric
            MetricDivergenceFreeProcess<2> divergencefree_process = MetricDivergenceFreeProcess<2>(r_model_part, this_parameters);
            divergencefree_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 4;
            ref_metric[1] = 4;
            ref_metric[2] = 0;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the divergencefree metric process with mean distribution strategy
        * Test triangle
        */

        KRATOS_TEST_CASE_IN_SUITE(DivergenceFreeMetricProcess2, ExaquteSandboxApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            Parameters this_parameters = Parameters(R"(
            {
                "minimal_size"                        : 0.01,
                "maximal_size"                        : 10.0,
                "refinement_strategy"                 : "mean_distribution_strategy",
                "reference_variable_name"             : "DISTANCE"
            })"
            );

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            Create2DGeometry(r_model_part, "Element2D3N");

            // Set DISTANCE variable
            for (std::size_t i_elem = 0; i_elem < r_model_part.Elements().size(); ++i_elem) {
                auto it_elem = r_model_part.Elements().begin() + i_elem;
                it_elem->SetValue(DISTANCE, 1.0);
            }

            // Set other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_2D, ZeroVector(3));
            }

            // Compute metric
            MetricDivergenceFreeProcess<2> divergencefree_process = MetricDivergenceFreeProcess<2>(r_model_part, this_parameters);
            divergencefree_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 1.23457;
            ref_metric[1] = 1.23457;
            ref_metric[2] = 0;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the divergencefree metric process with maximum strategy
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(DivergenceFreeMetricProcess3, ExaquteSandboxApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            Parameters this_parameters = Parameters(R"(
            {
                "minimal_size"                        : 0.01,
                "maximal_size"                        : 10.0,
                "refinement_strategy"                 : "maximum_strategy",
                "reference_variable_name"             : "DISTANCE"
            })"
            );

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            Create3DGeometry(r_model_part, "Element3D4N");

            // Set DISTANCE variable
            for (std::size_t i_elem = 0; i_elem < r_model_part.Elements().size(); ++i_elem) {
                auto it_elem = r_model_part.Elements().begin() + i_elem;
                it_elem->SetValue(DISTANCE, 1.0);
            }

            // Set other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_3D, ZeroVector(6));
            }

            // Compute metric
            MetricDivergenceFreeProcess<3> divergencefree_process = MetricDivergenceFreeProcess<3>(r_model_part, this_parameters);
            divergencefree_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 4;
            ref_metric[1] = 4;
            ref_metric[2] = 4;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the divergencefree metric process with mean distribution strategy
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(DivergenceFreeMetricProcess4, ExaquteSandboxApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            Parameters this_parameters = Parameters(R"(
            {
                "minimal_size"                        : 0.01,
                "maximal_size"                        : 10.0,
                "refinement_strategy"                 : "mean_distribution_strategy",
                "reference_variable_name"             : "DISTANCE"
            })"
            );

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            Create3DGeometry(r_model_part, "Element3D4N");

            // Set DISTANCE variable
            for (std::size_t i_elem = 0; i_elem < r_model_part.Elements().size(); ++i_elem) {
                auto it_elem = r_model_part.Elements().begin() + i_elem;
                it_elem->SetValue(DISTANCE, 1.0);
            }

            // Set other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_3D, ZeroVector(6));
            }

            // Compute metric
            MetricDivergenceFreeProcess<3> divergencefree_process = MetricDivergenceFreeProcess<3>(r_model_part, this_parameters);
            divergencefree_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 1.23457;
            ref_metric[1] = 1.23457;
            ref_metric[2] = 1.23457;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
        }

                /**
        * Checks the correct work of the divergencefree metric process with global tolerance strategy
        * Test triangle
        */

        KRATOS_TEST_CASE_IN_SUITE(DivergenceFreeMetricProcess5, ExaquteSandboxApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            Parameters this_parameters = Parameters(R"(
            {
                "minimal_size"                        : 0.01,
                "maximal_size"                        : 10.0,
                "refinement_strategy"                 : "global_tolerance_strategy",
                "reference_variable_name"             : "DISTANCE"
            })"
            );

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            Create2DGeometry(r_model_part, "Element2D3N");

            // Set DISTANCE variable
            for (std::size_t i_elem = 0; i_elem < r_model_part.Elements().size(); ++i_elem) {
                auto it_elem = r_model_part.Elements().begin() + i_elem;
                it_elem->SetValue(DISTANCE, 1.0);
            }

            // Set other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_2D, ZeroVector(3));
            }

            // Compute metric
            MetricDivergenceFreeProcess<2> divergencefree_process = MetricDivergenceFreeProcess<2>(r_model_part, this_parameters);
            divergencefree_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 10.0;
            ref_metric[1] = 10.0;
            ref_metric[2] = 0;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the divergencefree metric process with global tolerance strategy
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(DivergenceFreeMetricProcess6, ExaquteSandboxApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            Parameters this_parameters = Parameters(R"(
            {
                "minimal_size"                        : 0.01,
                "maximal_size"                        : 10.0,
                "refinement_strategy"                 : "global_tolerance_strategy",
                "reference_variable_name"             : "DISTANCE"
            })"
            );

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            Create3DGeometry(r_model_part, "Element3D4N");

            // Set DISTANCE variable
            for (std::size_t i_elem = 0; i_elem < r_model_part.Elements().size(); ++i_elem) {
                auto it_elem = r_model_part.Elements().begin() + i_elem;
                it_elem->SetValue(DISTANCE, 1.0);
            }

            // Set other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_3D, ZeroVector(6));
            }

            // Compute metric
            MetricDivergenceFreeProcess<3> divergencefree_process = MetricDivergenceFreeProcess<3>(r_model_part, this_parameters);
            divergencefree_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 10;
            ref_metric[1] = 10;
            ref_metric[2] = 10;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
        }

    } // namespace Testing
}  // namespace Kratos.
