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
#include "testing/testing.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "containers/model.h"
#include "meshing_application_variables.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/compute_nodal_gradient_process.h"
#include "custom_processes/metrics_hessian_process.h"
#include "custom_processes/metrics_levelset_process.h"
#include "custom_processes/metrics_error_process.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;

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

        /**
        * Checks the correct work of the level set metric process
        * Test triangle
        */

        KRATOS_TEST_CASE_IN_SUITE(LevelSetMetricProcess1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

            // Set DISTANCE and other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(NODAL_AREA, 0.0);
                it_node->SetValue(METRIC_TENSOR_2D, ZeroVector(3));
            }

            typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> GradientType;
            GradientType gradient_process = GradientType(r_model_part, DISTANCE, DISTANCE_GRADIENT, NODAL_AREA);
            gradient_process.Execute();

            // Compute metric
            ComputeLevelSetSolMetricProcess<2> level_set_process = ComputeLevelSetSolMetricProcess<2>(r_model_part);
            level_set_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 100;
            ref_metric[1] = 100;
            ref_metric[2] = 0;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the nodal gradient compute
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(LevelSetMetricProcess2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create3DGeometry(r_model_part, "Element3D4N");

            // Set DISTANCE and other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(NODAL_AREA, 0.0);
                it_node->SetValue(METRIC_TENSOR_3D, ZeroVector(6));
            }

            // Compute gradient
            typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> GradientType;
            GradientType gradient_process = GradientType(r_model_part, DISTANCE, DISTANCE_GRADIENT, NODAL_AREA);
            gradient_process.Execute();

            // Compute metric
            ComputeLevelSetSolMetricProcess<3> level_set_process = ComputeLevelSetSolMetricProcess<3>(r_model_part);
            level_set_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 100;
            ref_metric[1] = 100;
            ref_metric[2] = 100;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(3)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(9)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(10)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(11)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(12)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the hessian metric process
        * Test triangle
        */

        KRATOS_TEST_CASE_IN_SUITE(HessianMetricProcess1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N");

            // Set DISTANCE and other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_2D, ZeroVector(3));
            }

            // Compute metric
            Parameters parameters = Parameters(R"({"enforce_anisotropy_relative_variable" : true})");
            auto hessian_process = ComputeHessianSolMetricProcess(r_model_part, DISTANCE, parameters);
            hessian_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 100;
            ref_metric[1] = 100;
            ref_metric[2] = 0;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(6)->GetValue(METRIC_TENSOR_2D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the nodal gradient compute
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(HessianMetricProcess2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISTANCE);
            r_model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            CppTestsUtilities::Create3DGeometry(r_model_part, "Element3D4N");

            // Set DISTANCE and other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISTANCE) = (it_node->X() == 1.0) ? 0.0 : 1.0;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_3D, ZeroVector(6));
            }

            // Compute metric
            Parameters parameters = Parameters(R"({"enforce_anisotropy_relative_variable" : true})");
            auto hessian_process = ComputeHessianSolMetricProcess(r_model_part, DISTANCE, parameters);
            hessian_process.Execute();

//             // DEBUG
//             GiDIODebugMetric(r_model_part);

            const double tolerance = 1.0e-4;
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 100;
            ref_metric[1] = 100;
            ref_metric[2] = 100;
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(1)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(2)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(3)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(5)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(9)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(10)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(11)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
            KRATOS_CHECK_LESS_EQUAL(norm_2(r_model_part.pGetNode(12)->GetValue(METRIC_TENSOR_3D) - ref_metric), tolerance);
        }

        /**
        * Checks the correct work of the SPR metric process
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(SPRMetricProcess1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 2);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElasticPlaneStrain2DLaw"))
                return void();

            CppTestsUtilities::Create2DGeometry(r_model_part, "SmallDisplacementElement2D3N", false);

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            Properties::Pointer p_elem_prop = r_model_part.pGetProperties(0);

            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            p_elem_prop->SetValue(YOUNG_MODULUS, 1.0);
            p_elem_prop->SetValue(POISSON_RATIO, 0.0);

            for (auto& ielem : r_model_part.Elements()) {
                ielem.Initialize();
                ielem.InitializeSolutionStep(process_info);
            }

            // Set DISPLACEMENT_X and other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISPLACEMENT_X) = (it_node->X() == 1.0) ? 0.5 : 0.0;
                it_node->Coordinates()[0] += (it_node->X() == 1.0) ? 0.5 : 0.0;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_2D, ZeroVector(3));
            }

            // Compute error
            process_info[ERROR_OVERALL] = 0.122409;
            process_info[ENERGY_NORM_OVERALL] = 0.257196;
            for (auto& elem : r_model_part.Elements())
                elem.SetValue(ELEMENT_ERROR, 0.025);

            // Compute metric
            MetricErrorProcess<2> metric_process = MetricErrorProcess<2>(r_model_part);
            metric_process.Execute();

//             // DEBUG
//             GiDIODebugMetricSPR(r_model_part);

            const double tolerance = 1.0e-4;
            KRATOS_CHECK_LESS_EQUAL(r_model_part.pGetNode(2)->GetValue(METRIC_SCALAR) - std::sqrt(1.0/246.507)/r_model_part.pGetNode(2)->GetValue(METRIC_SCALAR), tolerance);
            KRATOS_CHECK_LESS_EQUAL(r_model_part.pGetNode(3)->GetValue(METRIC_SCALAR) - std::sqrt(1.0/246.507)/r_model_part.pGetNode(3)->GetValue(METRIC_SCALAR), tolerance);
        }

        /**
        * Checks the correct work of the nodal SPR compute
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(SPRMetricProcess2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info.SetValue(DOMAIN_SIZE, 3);
            process_info.SetValue(STEP, 1);
            process_info.SetValue(NL_ITERATION_NUMBER, 1);

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();

            CppTestsUtilities::Create3DGeometry(r_model_part, "SmallDisplacementElement3D4N", false);

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            Properties::Pointer p_elem_prop = r_model_part.pGetProperties(0);

            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            p_elem_prop->SetValue(YOUNG_MODULUS, 1.0);
            p_elem_prop->SetValue(POISSON_RATIO, 0.0);

            for (auto& ielem : r_model_part.Elements()) {
                ielem.Initialize();
                ielem.InitializeSolutionStep(process_info);
            }

            // Set DISPLACEMENT_X and other variables
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->FastGetSolutionStepValue(DISPLACEMENT_X) = (it_node->X() == 1.0) ? 0.5 : 0.0;
                it_node->Coordinates()[0] += (it_node->X() == 1.0) ? 0.5 : 0.0;
                it_node->SetValue(NODAL_H, 1.0);
                it_node->SetValue(METRIC_TENSOR_3D, ZeroVector(6));
            }

            // Compute error
            process_info[ERROR_OVERALL] = 0.0223607;
            process_info[ENERGY_NORM_OVERALL] = 0.148492;
            for (auto& elem : r_model_part.Elements())
                elem.SetValue(ELEMENT_ERROR, 0.025);

            // Compute metric
            MetricErrorProcess<3> metric_process = MetricErrorProcess<3>(r_model_part);
            metric_process.Execute();

//             // DEBUG
//             GiDIODebugMetricSPR(r_model_part);

            const double tolerance = 1.0e-4;
            KRATOS_CHECK_LESS_EQUAL(r_model_part.pGetNode(3)->GetValue(METRIC_SCALAR) - std::sqrt(1.0/(0.4807502774165066 * 4190.45))/r_model_part.pGetNode(3)->GetValue(METRIC_SCALAR), tolerance);
            KRATOS_CHECK_LESS_EQUAL(r_model_part.pGetNode(6)->GetValue(METRIC_SCALAR) - std::sqrt(1.0/4190.45)/r_model_part.pGetNode(6)->GetValue(METRIC_SCALAR), tolerance);
            KRATOS_CHECK_LESS_EQUAL(r_model_part.pGetNode(7)->GetValue(METRIC_SCALAR) - std::sqrt(1.0/4190.45)/r_model_part.pGetNode(7)->GetValue(METRIC_SCALAR), tolerance);
            KRATOS_CHECK_LESS_EQUAL(r_model_part.pGetNode(8)->GetValue(METRIC_SCALAR) - std::sqrt(1.0/4190.45)/r_model_part.pGetNode(8)->GetValue(METRIC_SCALAR), tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
