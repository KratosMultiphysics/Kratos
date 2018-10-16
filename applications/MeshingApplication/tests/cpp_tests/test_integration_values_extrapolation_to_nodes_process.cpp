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
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "includes/gid_io.h"
#include "meshing_application.h"

/* Processes */
#include "custom_processes/integration_values_extrapolation_to_nodes_process.h"
#include "includes/mat_variables.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

        void GiDIODebugInternalExtrapolation(ModelPart& ThisModelPart, const std::string name = "")
        {
            GidIO<> gid_io("TEST_INTERNAL_EXTRAPOLATION"+name, GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
            auto this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
            gid_io.PrintOnGaussPoints(this_var, ThisModelPart, label);
            gid_io.WriteNodalResultsNonHistorical(this_var, ThisModelPart.Nodes(), label);
            gid_io.WriteNodalResultsNonHistorical(NODAL_AREA, ThisModelPart.Nodes(), label);
        }
        
        void Create2DModelPartForExtrapolation(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_elem_prop = ThisModelPart.pGetProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;

            Element::Pointer p_elem_0 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 1, PointerVector<NodeType>{element_nodes_0}, p_elem_prop);
            Element::Pointer p_elem_1 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 2, PointerVector<NodeType>{element_nodes_1}, p_elem_prop);
            Element::Pointer p_elem_2 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 3, PointerVector<NodeType>{element_nodes_2}, p_elem_prop);
            Element::Pointer p_elem_3 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement2D3N", 4, PointerVector<NodeType>{element_nodes_3}, p_elem_prop);

            // Initialize Elements
            p_elem_0->Initialize();
            p_elem_1->Initialize();
            p_elem_2->Initialize();
            p_elem_3->Initialize();
        }

        void Create3DModelPartForExtrapolation(ModelPart& ThisModelPart)
        {
            Properties::Pointer p_elem_prop = ThisModelPart.pGetProperties(0);

            // First we create the nodes
            NodeType::Pointer p_node_1 = ThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_2 = ThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_3 = ThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_4 = ThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_5 = ThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = ThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

            NodeType::Pointer p_node_7 = ThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_8 = ThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_9 = ThisModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_10 = ThisModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_11 = ThisModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_12 = ThisModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (4);
            element_nodes_0[0] = p_node_12;
            element_nodes_0[1] = p_node_10;
            element_nodes_0[2] = p_node_8;
            element_nodes_0[3] = p_node_9;

            std::vector<NodeType::Pointer> element_nodes_1 (4);
            element_nodes_1[0] = p_node_4;
            element_nodes_1[1] = p_node_6;
            element_nodes_1[2] = p_node_9;
            element_nodes_1[3] = p_node_7;

            std::vector<NodeType::Pointer> element_nodes_2 (4);
            element_nodes_2[0] = p_node_11;
            element_nodes_2[1] = p_node_7;
            element_nodes_2[2] = p_node_9;
            element_nodes_2[3] = p_node_8;

            std::vector<NodeType::Pointer> element_nodes_3 (4);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_3;
            element_nodes_3[2] = p_node_8;
            element_nodes_3[3] = p_node_6;

            std::vector<NodeType::Pointer> element_nodes_4 (4);
            element_nodes_4[0] = p_node_4;
            element_nodes_4[1] = p_node_6;
            element_nodes_4[2] = p_node_7;
            element_nodes_4[3] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_5 (4);
            element_nodes_5[0] = p_node_2;
            element_nodes_5[1] = p_node_3;
            element_nodes_5[2] = p_node_5;
            element_nodes_5[3] = p_node_6;

            std::vector<NodeType::Pointer> element_nodes_6 (4);
            element_nodes_6[0] = p_node_10;
            element_nodes_6[1] = p_node_9;
            element_nodes_6[2] = p_node_6;
            element_nodes_6[3] = p_node_8;

            std::vector<NodeType::Pointer> element_nodes_7 (4);
            element_nodes_7[0] = p_node_7;
            element_nodes_7[1] = p_node_8;
            element_nodes_7[2] = p_node_3;
            element_nodes_7[3] = p_node_6;

            std::vector<NodeType::Pointer> element_nodes_8 (4);
            element_nodes_8[0] = p_node_7;
            element_nodes_8[1] = p_node_8;
            element_nodes_8[2] = p_node_6;
            element_nodes_8[3] = p_node_9;

            std::vector<NodeType::Pointer> element_nodes_9 (4);
            element_nodes_9[0] = p_node_4;
            element_nodes_9[1] = p_node_1;
            element_nodes_9[2] = p_node_6;
            element_nodes_9[3] = p_node_3;

            std::vector<NodeType::Pointer> element_nodes_10 (4);
            element_nodes_10[0] = p_node_9;
            element_nodes_10[1] = p_node_12;
            element_nodes_10[2] = p_node_11;
            element_nodes_10[3] = p_node_8;

            std::vector<NodeType::Pointer> element_nodes_11 (4);
            element_nodes_11[0] = p_node_3;
            element_nodes_11[1] = p_node_2;
            element_nodes_11[2] = p_node_1;
            element_nodes_11[3] = p_node_6;

            Element::Pointer p_elem_0 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 1, PointerVector<NodeType>{element_nodes_0}, p_elem_prop);
            Element::Pointer p_elem_1 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 2, PointerVector<NodeType>{element_nodes_1}, p_elem_prop);
            Element::Pointer p_elem_2 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 3, PointerVector<NodeType>{element_nodes_2}, p_elem_prop);
            Element::Pointer p_elem_3 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 4, PointerVector<NodeType>{element_nodes_3}, p_elem_prop);
            Element::Pointer p_elem_4 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 5, PointerVector<NodeType>{element_nodes_4}, p_elem_prop);
            Element::Pointer p_elem_5 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 6, PointerVector<NodeType>{element_nodes_5}, p_elem_prop);
            Element::Pointer p_elem_6 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 7, PointerVector<NodeType>{element_nodes_6}, p_elem_prop);
            Element::Pointer p_elem_7 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 8, PointerVector<NodeType>{element_nodes_7}, p_elem_prop);
            Element::Pointer p_elem_8 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 9, PointerVector<NodeType>{element_nodes_8}, p_elem_prop);
            Element::Pointer p_elem_9 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 10, PointerVector<NodeType>{element_nodes_9}, p_elem_prop);
            Element::Pointer p_elem_10 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 11, PointerVector<NodeType>{element_nodes_10}, p_elem_prop);
            Element::Pointer p_elem_11 = ThisModelPart.CreateNewElement("UpdatedLagrangianElement3D4N", 12, PointerVector<NodeType>{element_nodes_11}, p_elem_prop);

            // Initialize Elements
            p_elem_0->Initialize();
            p_elem_1->Initialize();
            p_elem_2->Initialize();
            p_elem_3->Initialize();
            p_elem_4->Initialize();
            p_elem_5->Initialize();
            p_elem_6->Initialize();
            p_elem_7->Initialize();
            p_elem_8->Initialize();
            p_elem_9->Initialize();
            p_elem_10->Initialize();
            p_elem_11->Initialize();
        }

        /**
        * Checks the correct work of the internal variable extrapolation process
        * Test triangle 
        */

        KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessTriangle, KratosMeshingApplicationFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;
            
            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElasticPlaneStrain2DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            Create2DModelPartForExtrapolation(this_model_part);

            // Compute extrapolation
            Parameters extrapolation_parameters = Parameters(R"(
            {
                "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            })");
            IntegrationValuesExtrapolationToNodesProcess extrapolation_process = IntegrationValuesExtrapolationToNodesProcess(this_model_part, extrapolation_parameters);
            extrapolation_process.ExecuteBeforeSolutionLoop();
            extrapolation_process.ExecuteFinalizeSolutionStep();

//             // DEBUG
//             GiDIODebugInternalExtrapolation(this_model_part, "1");
            
            const double tolerance = 1.0e-8;
            auto this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
            for (auto& node : this_model_part.Nodes()) {
                KRATOS_CHECK_LESS_EQUAL(std::abs(node.GetValue(this_var) - 1.0), tolerance);
            }

            extrapolation_process.ExecuteFinalize();
        }
        
        /**
        * Checks the correct work of the internal variable extrapolation process
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(TestIntegrationValuesExtrapolationToNodesProcessTetra, KratosMeshingApplicationFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            Create3DModelPartForExtrapolation(this_model_part);

            // Compute extrapolation
            Parameters extrapolation_parameters = Parameters(R"(
            {
                "list_of_variables"          : ["REFERENCE_DEFORMATION_GRADIENT_DETERMINANT"]
            })");
            IntegrationValuesExtrapolationToNodesProcess extrapolation_process = IntegrationValuesExtrapolationToNodesProcess(this_model_part, extrapolation_parameters);
            extrapolation_process.ExecuteBeforeSolutionLoop();
            extrapolation_process.ExecuteFinalizeSolutionStep();

//             // DEBUG
//             GiDIODebugInternalExtrapolation(this_model_part, "2");

            const double tolerance = 1.0e-8;
            auto this_var = KratosComponents<Variable<double>>::Get("REFERENCE_DEFORMATION_GRADIENT_DETERMINANT");
            for (auto& node : this_model_part.Nodes()) {
                KRATOS_CHECK_LESS_EQUAL(std::abs(node.GetValue(this_var) - 1.0), tolerance);
            }

            extrapolation_process.ExecuteFinalize();
        }
    } // namespace Testing
}  // namespace Kratos.
#endif
