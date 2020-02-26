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
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/kratos_flags.h"
// #include "includes/gid_io.h"

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

            // First we create the nodes
            rModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            rModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            rModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Now we create the elements
            rModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 1, {{1,2,3}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 2, {{1,3,4}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 3, {{2,5,3}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 4, {{5,6,3}}, p_elem_prop);

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
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 1, {{12,10,8,9}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 2, {{4,6,9,7}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 3, {{11,7,9,8}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 4, {{5,3,8,6}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 5, {{4,6,7,3}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 6, {{2,3,5,6}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 7, {{10,9,6,8}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 8, {{7,8,3,6}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 9, {{7,8,6,9}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 10, {{4,1,6,3}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 11, {{9,12,11,8}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 12, {{3,2,1,6}}, p_elem_prop);

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

            // First we create the nodes
            rModelPart.CreateNewNode(1, 1.0000000000, 0.0000000000, 1.0000000000);
            rModelPart.CreateNewNode(2, 0.5000000000, 0.0000000000, 1.0000000000);
            rModelPart.CreateNewNode(3, 1.0000000000, 0.5000000000, 1.0000000000);
            rModelPart.CreateNewNode(4, 1.0000000000, 0.0000000000, 0.5000000000);
            rModelPart.CreateNewNode(5, 0.5000000000, 0.0000000000, 0.5000000000);
            rModelPart.CreateNewNode(6, 1.0000000000, 0.5000000000, 0.5000000000);
            rModelPart.CreateNewNode(7, 0.5000000000, 0.5000000000, 1.0000000000);
            rModelPart.CreateNewNode(8, 0.5000000000, 0.5000000000, 0.5000000000);
            rModelPart.CreateNewNode(9, 0.2019246055, 0.3959160307, 0.6930668948);
            rModelPart.CreateNewNode(10, 0.7019246055, 0.8959160307, 0.6930668948);
            rModelPart.CreateNewNode(11, 1.0000000000, 0.0000000000, 0.0000000000);
            rModelPart.CreateNewNode(12, 0.0000000000, 0.0000000000, 1.0000000000);
            rModelPart.CreateNewNode(13, 1.0000000000, 1.0000000000, 1.0000000000);
            rModelPart.CreateNewNode(14, 0.5000000000, 0.0000000000, 0.0000000000);
            rModelPart.CreateNewNode(15, 1.0000000000, 0.5000000000, 0.0000000000);
            rModelPart.CreateNewNode(16, 0.5000000000, 1.0000000000, 1.0000000000);
            rModelPart.CreateNewNode(17, 0.0000000000, 0.5000000000, 1.0000000000);
            rModelPart.CreateNewNode(18, 0.0000000000, 0.0000000000, 0.5000000000);
            rModelPart.CreateNewNode(19, 1.0000000000, 1.0000000000, 0.5000000000);
            rModelPart.CreateNewNode(20, 0.4038492111, 0.7918320615, 0.3861337896);
            rModelPart.CreateNewNode(21, 0.2019246055, 0.3959160307, 0.1930668948);
            rModelPart.CreateNewNode(22, 0.5000000000, 0.5000000000, 0.0000000000);
            rModelPart.CreateNewNode(23, 0.5000000000, 1.0000000000, 0.5000000000);
            rModelPart.CreateNewNode(24, 0.0000000000, 0.5000000000, 0.5000000000);
            rModelPart.CreateNewNode(25, 0.2019246055, 0.8959160307, 0.6930668948);
            rModelPart.CreateNewNode(26, 0.7019246055, 0.8959160307, 0.1930668948);
            rModelPart.CreateNewNode(27, 0.0000000000, 0.0000000000, 0.0000000000);
            rModelPart.CreateNewNode(28, 1.0000000000, 1.0000000000, 0.0000000000);
            rModelPart.CreateNewNode(29, 0.0000000000, 1.0000000000, 1.0000000000);
            rModelPart.CreateNewNode(30, 0.2019246055, 0.8959160307, 0.1930668948);
            rModelPart.CreateNewNode(31, 0.5000000000, 1.0000000000, 0.0000000000);
            rModelPart.CreateNewNode(32, 0.0000000000, 0.5000000000, 0.0000000000);
            rModelPart.CreateNewNode(33, 0.0000000000, 1.0000000000, 0.5000000000);
            rModelPart.CreateNewNode(34, 0.0000000000, 1.0000000000, 0.0000000000);

            // Now we create the elements
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 1,  {{27, 11, 28, 13, 14, 15, 22,  8,  6, 19}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 2,  {{34, 12, 27, 20, 24, 18, 32, 30,  9, 21}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 3,  {{27, 34, 20, 28, 32, 30, 21, 22, 31, 26}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 4,  {{34, 20, 28, 29, 30, 26, 31, 33, 25, 23}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 5,  {{34, 20, 29, 12, 30, 25, 33, 24,  9, 17}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 6,  {{20, 27, 28, 13, 21, 22, 26, 10,  8, 19}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 7,  {{28, 20, 13, 29, 26, 10, 19, 23, 25, 16}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 8,  {{20, 13, 29, 12, 10, 16, 25,  9,  7, 17}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 9,  {{27,  1, 11, 13,  5,  4, 14,  8,  3,  6}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 10, {{12, 13,  1, 27,  7,  3,  2, 18,  8,  5}}, p_elem_prop);
            rModelPart.CreateNewElement("UpdatedLagrangianElement3D10N", 11, {{12, 27, 20, 13, 18, 21,  9,  7,  8, 10}}, p_elem_prop);

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
