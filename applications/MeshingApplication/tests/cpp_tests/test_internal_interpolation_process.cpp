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
#include "containers/model.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "meshing_application_variables.h"

/* Processes */
#include "processes/find_nodal_h_process.h"
#include "custom_processes/mmg_process.h"
#include "utilities/cpp_tests_utilities.h"
#include "includes/mat_variables.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

//         void GiDIODebugInternalInterpolation(ModelPart& rModelPart, const std::string name = "")
//         {
//             GidIO<> gid_io("TEST_INTERNAL_INTERPOLATION_MMG"+name, GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             const Variable<double>& plastic_strain_variable = KratosComponents<Variable<double>>::Get("ACCUMULATED_PLASTIC_STRAIN");
//             gid_io.PrintOnGaussPoints(plastic_strain_variable, rModelPart, label);
//         }
//
//         void GiDIODebugInternalInterpolationElement(ModelPart& rModelPart, const std::string name = "")
//         {
//             GidIO<> gid_io("TEST_INTERNAL_INTERPOLATION_MMG"+name, GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
//             gid_io.PrintOnGaussPoints(this_var, rModelPart, label);
//         }

        /**
        * Checks the correct work of the internal variable interpolation process after remesh CPT test
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessCPT1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("SmallStrainJ2PlasticityPlaneStrain2DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("SmallStrainJ2PlasticityPlaneStrain2DLaw");
            const Variable<double>& plastic_strain_variable = KratosComponents<Variable<double>>::Get("ACCUMULATED_PLASTIC_STRAIN");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(r_model_part, "UpdatedLagrangianElement2D3N");

            // Set DISTANCE and other variables
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 0.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_2D, ref_metric);
            }

            // Set ACCUMULATED_PLASTIC_STRAIN on the GP
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(plastic_strain_variable, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "CPT",
                "internal_variable_interpolation_list" : ["ACCUMULATED_PLASTIC_STRAIN"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "pre1");

            MmgProcess<MMGLibrary::MMG2D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "1");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    double aux;
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(plastic_strain_variable, aux) - 1.0), tolerance);
                }
            }
        }
        /**
        * Checks the correct work of the internal variable interpolation process after remesh LST test
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessLST1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("SmallStrainJ2PlasticityPlaneStrain2DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("SmallStrainJ2PlasticityPlaneStrain2DLaw");
            const Variable<double>& plastic_strain_variable = KratosComponents<Variable<double>>::Get("ACCUMULATED_PLASTIC_STRAIN");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(r_model_part, "UpdatedLagrangianElement2D3N");

            // Set DISTANCE and other variables
            array_1d<double, 3> ref_metric(3);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 0.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_2D, ref_metric);
            }

            // Set ACCUMULATED_PLASTIC_STRAIN on the GP
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(plastic_strain_variable, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" : ["ACCUMULATED_PLASTIC_STRAIN"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "pre1");

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
            process.Execute();

            MmgProcess<MMGLibrary::MMG2D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "1");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    double aux;
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(plastic_strain_variable, aux) - 1.0), tolerance);
                }
            }
        }

        /**
        * Checks the correct work of the internal variable interpolation process after remesh CPT test
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessCPT2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("SmallStrainJ2Plasticity3DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("SmallStrainJ2Plasticity3DLaw");
            const Variable<double>& plastic_strain_variable = KratosComponents<Variable<double>>::Get("ACCUMULATED_PLASTIC_STRAIN");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create3DGeometry(r_model_part, "UpdatedLagrangianElement3D4N");

            // Set DISTANCE and other variables
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_3D, ref_metric);
            }

            // Set ACCUMULATED_PLASTIC_STRAIN on the GP
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(plastic_strain_variable, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "CPT",
                "internal_variable_interpolation_list" : ["ACCUMULATED_PLASTIC_STRAIN"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "pre2");

            MmgProcess<MMGLibrary::MMG3D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "2");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    double aux;
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(plastic_strain_variable, aux) - 1.0), tolerance);
                }
            }
        }

        /**
        * Checks the correct work of the internal variable interpolation process after remesh LST test
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessLST2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("SmallStrainJ2Plasticity3DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("SmallStrainJ2Plasticity3DLaw");
            const Variable<double>& plastic_strain_variable = KratosComponents<Variable<double>>::Get("ACCUMULATED_PLASTIC_STRAIN");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create3DGeometry(r_model_part, "UpdatedLagrangianElement3D4N");

            // Set DISTANCE and other variables
            array_1d<double, 6> ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_3D, ref_metric);
            }

            // Set ACCUMULATED_PLASTIC_STRAIN on the GP
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(plastic_strain_variable, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" : ["ACCUMULATED_PLASTIC_STRAIN"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "pre2");

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
            process.Execute();

            MmgProcess<MMGLibrary::MMG3D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(r_model_part, "2");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    double aux;
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(plastic_strain_variable, aux) - 1.0), tolerance);
                }
            }
        }

        /**
        * Checks the correct work of the internal variable interpolation process after remesh CPT test
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessElementsCPT1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElasticPlaneStrain2DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(r_model_part, "UpdatedLagrangianElement2D3N");

            // Set DISTANCE and other variables
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 0.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_2D, ref_metric);
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "CPT",
                "internal_variable_interpolation_list" : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "pre1");

            MmgProcess<MMGLibrary::MMG2D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "1");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the this_var
                std::vector<double> detF0_vector(integration_points_number);
                const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
                elem.GetValueOnIntegrationPoints(this_var,detF0_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++)
                    KRATOS_CHECK_LESS_EQUAL(std::abs(detF0_vector[i] - 1.0), tolerance);
            }
        }
        /**
        * Checks the correct work of the internal variable interpolation process after remesh LST test
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessElementsLST1, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElasticPlaneStrain2DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create2DGeometry(r_model_part, "UpdatedLagrangianElement2D3N");

            // Set DISTANCE and other variables
            array_1d<double, 3> ref_metric;
            ref_metric[0] = 1.0;
            ref_metric[1] = 1.0;
            ref_metric[2] = 0.0;
            for (std::size_t i_node = 0; i_node < r_model_part.Nodes().size(); ++i_node) {
                auto it_node = r_model_part.Nodes().begin() + i_node;
                it_node->SetValue(METRIC_TENSOR_2D, ref_metric);
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "pre1");

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
            process.Execute();

            MmgProcess<MMGLibrary::MMG2D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "1");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the this_var
                std::vector<double> detF0_vector(integration_points_number);
                const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
                elem.GetValueOnIntegrationPoints(this_var,detF0_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++)
                    KRATOS_CHECK_LESS_EQUAL(std::abs(detF0_vector[i] - 1.0), tolerance);
            }
        }

        /**
        * Checks the correct work of the internal variable interpolation process after remesh CPT test
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessElementsCPT2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create3DGeometry(r_model_part, "UpdatedLagrangianElement3D4N");

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
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "CPT",
                "internal_variable_interpolation_list" : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "pre2");

            MmgProcess<MMGLibrary::MMG3D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "2");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the this_var
                std::vector<double> detF0_vector(integration_points_number);
                const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
                elem.GetValueOnIntegrationPoints(this_var,detF0_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++)
                    KRATOS_CHECK_LESS_EQUAL(std::abs(detF0_vector[i] - 1.0), tolerance);
            }
        }

        /**
        * Checks the correct work of the internal variable interpolation process after remesh LST test
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(InternalInterpolationProcessElementsLST2, KratosMeshingApplicationFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& current_process_info = r_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = r_model_part.CreateNewProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            CppTestsUtilities::Create3DGeometry(r_model_part, "UpdatedLagrangianElement3D4N");

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
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "pre2");

            // Compute NodalH
            auto process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
            process.Execute();

            MmgProcess<MMGLibrary::MMG3D> mmg_process(r_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolationElement(r_model_part, "2");

            const double tolerance = 1.0e-4;
            for (auto& elem : r_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the this_var
                std::vector<double> detF0_vector(integration_points_number);
                const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
                elem.GetValueOnIntegrationPoints(this_var,detF0_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++)
                    KRATOS_CHECK_LESS_EQUAL(std::abs(detF0_vector[i] - 1.0), tolerance);
            }
        }
    } // namespace Testing
}  // namespace Kratos.
#endif
