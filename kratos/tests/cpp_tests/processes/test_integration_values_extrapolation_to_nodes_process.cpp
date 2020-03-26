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
#include "containers/model.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/integration_values_extrapolation_to_nodes_process.h"
#include "includes/mat_variables.h"

namespace Kratos
{
    namespace Testing
    {
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

//         void GiDIODebugInternalExtrapolation(ModelPart& rModelPart, const std::string name = "")
//         {
//             GidIO<> gid_io("TEST_INTERNAL_EXTRAPOLATION"+name, GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//             const int nl_iter = rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//             const double label = static_cast<double>(nl_iter);
//
//             gid_io.InitializeMesh(label);
//             gid_io.WriteMesh(rModelPart.GetMesh());
//             gid_io.FinalizeMesh();
//             gid_io.InitializeResults(label, rModelPart.GetMesh());
//             auto this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
//             gid_io.PrintOnGaussPoints(this_var, rModelPart, label);
//             gid_io.WriteNodalResultsNonHistorical(this_var, rModelPart.Nodes(), label);
//             gid_io.WriteNodalResultsNonHistorical(NODAL_AREA, rModelPart.Nodes(), label);
//         }

        void Create2DModelPartForExtrapolation(ModelPart& rModelPart)
        {
            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            CppTestsUtilities::Create2DGeometry(rModelPart, "UpdatedLagrangianElement2D3N", false);

            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            // Initialize Elements
            for (auto& r_elem : rModelPart.Elements())
                r_elem.Initialize();
        }

        void Create3DModelPartForExtrapolation(ModelPart& rModelPart)
        {
            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            CppTestsUtilities::Create3DGeometry(rModelPart, "UpdatedLagrangianElement3D4N", false);

            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            // Initialize Elements
            for (auto& r_elem : rModelPart.Elements())
                r_elem.Initialize();
        }

        void CreateQuadratic3DModelPartForExtrapolation(ModelPart& rModelPart)
        {
            Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

            CppTestsUtilities::Create3DQuadraticGeometry(rModelPart, "UpdatedLagrangianElement3D10N", false);

            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            // Initialize Elements
            for (auto& r_elem : rModelPart.Elements())
                r_elem.Initialize();
        }

        /**
        * Checks the correct work of the internal variable extrapolation process
        * Test triangle
        */
        KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessTriangle, KratosCoreFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 2;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElasticPlaneStrain2DLaw"))
                return void();

            Create2DModelPartForExtrapolation(r_model_part);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // Compute extrapolation
            Parameters extrapolation_parameters = Parameters(R"(
            {
                "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            })");
            IntegrationValuesExtrapolationToNodesProcess extrapolation_process(r_model_part, extrapolation_parameters);
            extrapolation_process.ExecuteBeforeSolutionLoop();
            extrapolation_process.ExecuteFinalizeSolutionStep();

//             // DEBUG
//             GiDIODebugInternalExtrapolation(r_model_part, "1");

            const double tolerance = 1.0e-8;
            const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
            for (auto& node : r_model_part.Nodes()) {
                KRATOS_CHECK_LESS_EQUAL(std::abs(node.GetValue(this_var) - 1.0), tolerance);
            }

            extrapolation_process.ExecuteFinalize();
        }

        /**
        * Checks the correct work of the internal variable extrapolation process
        * Test tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessTetra, KratosCoreFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 3;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();

            Create3DModelPartForExtrapolation(r_model_part);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // Compute extrapolation
            Parameters extrapolation_parameters = Parameters(R"(
            {
                "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            })");
            IntegrationValuesExtrapolationToNodesProcess extrapolation_process(r_model_part, extrapolation_parameters);
            extrapolation_process.ExecuteBeforeSolutionLoop();
            extrapolation_process.ExecuteFinalizeSolutionStep();

//             // DEBUG
//             GiDIODebugInternalExtrapolation(r_model_part, "2");

            const double tolerance = 1.0e-8;
            const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
            for (auto& node : r_model_part.Nodes()) {
                KRATOS_CHECK_LESS_EQUAL(std::abs(node.GetValue(this_var) - 1.0), tolerance);
            }

            extrapolation_process.ExecuteFinalize();
        }

        /**
        * Checks the correct work of the internal variable extrapolation process
        * Test quadratic tetrahedra
        */
        KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessQuadTetra, KratosCoreFastSuite)
        {
            Model this_model;
            ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            r_current_process_info[DOMAIN_SIZE] = 3;

            r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();

            CreateQuadratic3DModelPartForExtrapolation(r_model_part);

            auto& process_info = r_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // Compute extrapolation
            Parameters extrapolation_parameters = Parameters(R"(
            {
                "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            })");
            IntegrationValuesExtrapolationToNodesProcess extrapolation_process(r_model_part, extrapolation_parameters);
            extrapolation_process.ExecuteBeforeSolutionLoop();
            extrapolation_process.ExecuteFinalizeSolutionStep();

//             // DEBUG
//             GiDIODebugInternalExtrapolation(r_model_part, "3");

            const double tolerance = 1.0e-6;
            const auto& this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
            for (auto& node : r_model_part.Nodes()) {
                KRATOS_CHECK_LESS_EQUAL(std::abs(node.GetValue(this_var) - 1.0), tolerance);
            }

            extrapolation_process.ExecuteFinalize();
        }
    } // namespace Testing
}  // namespace Kratos.
