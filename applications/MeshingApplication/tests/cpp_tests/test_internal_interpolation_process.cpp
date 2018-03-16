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

// Elements and CL
#include "custom_elements/small_displacement_bbar.h"
#include "custom_constitutive/linear_j2_plasticity_plane_strain_2d.h"
#include "custom_constitutive/linear_j2_plasticity_3d.h"

/* Processes */
#include "processes/find_nodal_h_process.h"
#include "custom_processes/mmg_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        
        void GiDIODebugInternalInterpolation(ModelPart& ThisModelPart, const std::string name = "")
        {
            GidIO<> gid_io("TEST_INTERNAL_INTERPOLATION_MMG"+name, GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
            const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
            const double label = static_cast<double>(nl_iter);

            gid_io.InitializeMesh(label);
            gid_io.WriteMesh(ThisModelPart.GetMesh());
            gid_io.FinalizeMesh();
            gid_io.InitializeResults(label, ThisModelPart.GetMesh());
            gid_io.PrintOnGaussPoints(PLASTIC_STRAIN, ThisModelPart, label);
        }
        
        /**
        * Checks the correct work of the internal variable interpolation process after remesh CPT test
        * Test triangle
        */

        KRATOS_TEST_CASE_IN_SUITE(TestInternalInterpolationProcessCPT1, KratosInternalInterpolationProcessFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;

            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            auto p_this_law = Kratos::make_shared<LinearJ2PlasticityPlaneStrain2D>();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // First we create the nodes
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_0( element_nodes_0 );

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle2D3 <NodeType> triangle_1( element_nodes_1 );

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_2( element_nodes_2 );

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_3( element_nodes_3 );

            // Creating pointers
            auto p_triangle_0 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_0);
            auto p_triangle_1 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_1);
            auto p_triangle_2 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_2);
            auto p_triangle_3 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_3);

            Element::Pointer p_elem_0 = Kratos::make_shared<SmallDisplacementBbar>(1, p_triangle_0, p_elem_prop);;
            Element::Pointer p_elem_1 = Kratos::make_shared<SmallDisplacementBbar>(2, p_triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = Kratos::make_shared<SmallDisplacementBbar>(3, p_triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = Kratos::make_shared<SmallDisplacementBbar>(4, p_triangle_3, p_elem_prop);

            // Adding to the model part
            this_model_part.AddElement(p_elem_0);
            this_model_part.AddElement(p_elem_1);
            this_model_part.AddElement(p_elem_2);
            this_model_part.AddElement(p_elem_3);

            // Initialize Elements
            p_elem_0->Initialize();
            p_elem_1->Initialize();
            p_elem_2->Initialize();
            p_elem_3->Initialize();

            // Set DISTANCE and other variables
            Vector ref_metric(3);
            ref_metric[0] = 1.0;
            ref_metric[1] = 0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->SetValue(MMG_METRIC, ref_metric);
            }

            // Set PLASTIC_STRAIN on the GP
            for (auto& elem : this_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(PLASTIC_STRAIN, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "CPT",
                "internal_variable_interpolation_list" : ["PLASTIC_STRAIN"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "pre1");

            MmgProcess<2> mmg_process = MmgProcess<2>(this_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "1");

            const double tolerance = 1.0e-4;
            for (auto& elem : this_model_part.Elements()) {
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
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(PLASTIC_STRAIN, aux) - 1.0), tolerance);
                }
            }
        }
        /**
        * Checks the correct work of the internal variable interpolation process after remesh LST test
        * Test triangle 
        */

        KRATOS_TEST_CASE_IN_SUITE(TestInternalInterpolationProcessLST1, KratosInternalInterpolationProcessFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 2;
            
            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);
            
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            auto p_this_law = Kratos::make_shared<LinearJ2PlasticityPlaneStrain2D>();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5, 2.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6, 2.0 , 1.0 , 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (3);
            element_nodes_0[0] = p_node_1;
            element_nodes_0[1] = p_node_2;
            element_nodes_0[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_0( element_nodes_0 );
            
            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle2D3 <NodeType> triangle_1( element_nodes_1 );
            
            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_2( element_nodes_2 );
            
            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_3( element_nodes_3 );

            // Creating pointers
            auto p_triangle_0 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_0);
            auto p_triangle_1 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_1);
            auto p_triangle_2 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_2);
            auto p_triangle_3 = Kratos::make_shared<Triangle2D3 <NodeType>> (triangle_3);

            Element::Pointer p_elem_0 = Kratos::make_shared<SmallDisplacementBbar>(1, p_triangle_0, p_elem_prop);;
            Element::Pointer p_elem_1 = Kratos::make_shared<SmallDisplacementBbar>(2, p_triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = Kratos::make_shared<SmallDisplacementBbar>(3, p_triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = Kratos::make_shared<SmallDisplacementBbar>(4, p_triangle_3, p_elem_prop);
            
            // Adding to the model part
            this_model_part.AddElement(p_elem_0);
            this_model_part.AddElement(p_elem_1);
            this_model_part.AddElement(p_elem_2);
            this_model_part.AddElement(p_elem_3);

            // Initialize Elements
            p_elem_0->Initialize();
            p_elem_1->Initialize();
            p_elem_2->Initialize();
            p_elem_3->Initialize();
            
            // Set DISTANCE and other variables
            Vector ref_metric(3);
            ref_metric[0] = 1.0;
            ref_metric[1] = 0;
            ref_metric[2] = 1.0;
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->SetValue(MMG_METRIC, ref_metric);
            }
                         
            // Set PLASTIC_STRAIN on the GP
            for (auto& elem : this_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(PLASTIC_STRAIN, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
            "framework"                            : "Lagrangian",
            "internal_variables_parameters"        :
            {
                "interpolation_type"                   : "LST",
                "internal_variable_interpolation_list" : ["PLASTIC_STRAIN"]
            },
            "echo_level" : 0
            })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "pre1");

            // Compute NodalH
            FindNodalHProcess process = FindNodalHProcess(this_model_part);
            process.Execute();

            MmgProcess<2> mmg_process = MmgProcess<2>(this_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "1");
            
            const double tolerance = 1.0e-4;
            for (auto& elem : this_model_part.Elements()) {
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
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(PLASTIC_STRAIN, aux) - 1.0), tolerance);
                }
            }
        }
        
        /**
        * Checks the correct work of the internal variable interpolation process after remesh CPT test
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(TestInternalInterpolationProcessCPT2, KratosInternalInterpolationProcessFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            auto p_this_law = Kratos::make_shared<LinearJ2Plasticity3D>();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;

            // First we create the nodes
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_9 = this_model_part.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_10 = this_model_part.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_11 = this_model_part.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_12 = this_model_part.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (4);
            element_nodes_0[0] = p_node_12;
            element_nodes_0[1] = p_node_10;
            element_nodes_0[2] = p_node_8;
            element_nodes_0[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_0( element_nodes_0 );

            std::vector<NodeType::Pointer> element_nodes_1 (4);
            element_nodes_1[0] = p_node_4;
            element_nodes_1[1] = p_node_6;
            element_nodes_1[2] = p_node_9;
            element_nodes_1[3] = p_node_7;
            Tetrahedra3D4 <NodeType> tetrahedra_1( element_nodes_1 );

            std::vector<NodeType::Pointer> element_nodes_2 (4);
            element_nodes_2[0] = p_node_11;
            element_nodes_2[1] = p_node_7;
            element_nodes_2[2] = p_node_9;
            element_nodes_2[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_2( element_nodes_2 );

            std::vector<NodeType::Pointer> element_nodes_3 (4);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_3;
            element_nodes_3[2] = p_node_8;
            element_nodes_3[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_3( element_nodes_3 );

            std::vector<NodeType::Pointer> element_nodes_4 (4);
            element_nodes_4[0] = p_node_4;
            element_nodes_4[1] = p_node_6;
            element_nodes_4[2] = p_node_7;
            element_nodes_4[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_4( element_nodes_4 );

            std::vector<NodeType::Pointer> element_nodes_5 (4);
            element_nodes_5[0] = p_node_2;
            element_nodes_5[1] = p_node_3;
            element_nodes_5[2] = p_node_5;
            element_nodes_5[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_5( element_nodes_5 );

            std::vector<NodeType::Pointer> element_nodes_6 (4);
            element_nodes_6[0] = p_node_10;
            element_nodes_6[1] = p_node_9;
            element_nodes_6[2] = p_node_6;
            element_nodes_6[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_6( element_nodes_6 );

            std::vector<NodeType::Pointer> element_nodes_7 (4);
            element_nodes_7[0] = p_node_7;
            element_nodes_7[1] = p_node_8;
            element_nodes_7[2] = p_node_3;
            element_nodes_7[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_7( element_nodes_7 );

            std::vector<NodeType::Pointer> element_nodes_8 (4);
            element_nodes_8[0] = p_node_7;
            element_nodes_8[1] = p_node_8;
            element_nodes_8[2] = p_node_6;
            element_nodes_8[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_8( element_nodes_8 );

            std::vector<NodeType::Pointer> element_nodes_9 (4);
            element_nodes_9[0] = p_node_4;
            element_nodes_9[1] = p_node_1;
            element_nodes_9[2] = p_node_6;
            element_nodes_9[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_9( element_nodes_9 );

            std::vector<NodeType::Pointer> element_nodes_10 (4);
            element_nodes_10[0] = p_node_9;
            element_nodes_10[1] = p_node_12;
            element_nodes_10[2] = p_node_11;
            element_nodes_10[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_10( element_nodes_10 );

            std::vector<NodeType::Pointer> element_nodes_11 (4);
            element_nodes_11[0] = p_node_3;
            element_nodes_11[1] = p_node_2;
            element_nodes_11[2] = p_node_1;
            element_nodes_11[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_11( element_nodes_11 );

            // Creating pointers
            auto p_tetrahedra_0 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_0);
            auto p_tetrahedra_1 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_1);
            auto p_tetrahedra_2 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_2);
            auto p_tetrahedra_3 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_3);
            auto p_tetrahedra_4 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_4);
            auto p_tetrahedra_5 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_5);
            auto p_tetrahedra_6 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_6);
            auto p_tetrahedra_7 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_7);
            auto p_tetrahedra_8 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_8);
            auto p_tetrahedra_9 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_9);
            auto p_tetrahedra_10 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_10);
            auto p_tetrahedra_11 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_11);

            Element::Pointer p_elem_0 = Kratos::make_shared<SmallDisplacementBbar>(1, p_tetrahedra_0, p_elem_prop);
            Element::Pointer p_elem_1 = Kratos::make_shared<SmallDisplacementBbar>(2, p_tetrahedra_1, p_elem_prop);
            Element::Pointer p_elem_2 = Kratos::make_shared<SmallDisplacementBbar>(3, p_tetrahedra_2, p_elem_prop);
            Element::Pointer p_elem_3 = Kratos::make_shared<SmallDisplacementBbar>(4, p_tetrahedra_3, p_elem_prop);
            Element::Pointer p_elem_4 = Kratos::make_shared<SmallDisplacementBbar>(5, p_tetrahedra_4, p_elem_prop);
            Element::Pointer p_elem_5 = Kratos::make_shared<SmallDisplacementBbar>(6, p_tetrahedra_5, p_elem_prop);
            Element::Pointer p_elem_6 = Kratos::make_shared<SmallDisplacementBbar>(7, p_tetrahedra_6, p_elem_prop);
            Element::Pointer p_elem_7 = Kratos::make_shared<SmallDisplacementBbar>(8, p_tetrahedra_7, p_elem_prop);
            Element::Pointer p_elem_8 = Kratos::make_shared<SmallDisplacementBbar>(9, p_tetrahedra_8, p_elem_prop);
            Element::Pointer p_elem_9 = Kratos::make_shared<SmallDisplacementBbar>(10, p_tetrahedra_9, p_elem_prop);
            Element::Pointer p_elem_10 = Kratos::make_shared<SmallDisplacementBbar>(11, p_tetrahedra_10, p_elem_prop);
            Element::Pointer p_elem_11 = Kratos::make_shared<SmallDisplacementBbar>(12, p_tetrahedra_11, p_elem_prop);

            // Adding to the model part
            this_model_part.AddElement(p_elem_0);
            this_model_part.AddElement(p_elem_1);
            this_model_part.AddElement(p_elem_2);
            this_model_part.AddElement(p_elem_3);
            this_model_part.AddElement(p_elem_4);
            this_model_part.AddElement(p_elem_5);
            this_model_part.AddElement(p_elem_6);
            this_model_part.AddElement(p_elem_7);
            this_model_part.AddElement(p_elem_8);
            this_model_part.AddElement(p_elem_9);
            this_model_part.AddElement(p_elem_10);
            this_model_part.AddElement(p_elem_11);

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

            // Set DISTANCE and other variables
            Vector ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[3] = 1.0;
            ref_metric[5] = 1.0;
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->SetValue(MMG_METRIC, ref_metric);
            }

            // Set PLASTIC_STRAIN on the GP
            for (auto& elem : this_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(PLASTIC_STRAIN, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
                "framework"                            : "Lagrangian",
                "internal_variables_parameters"        :
                {
                    "interpolation_type"                   : "CPT",
                    "internal_variable_interpolation_list" : ["PLASTIC_STRAIN"]
                },
                "echo_level" : 0
                })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "pre2");

            MmgProcess<3> mmg_process = MmgProcess<3>(this_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "2");

            const double tolerance = 1.0e-4;
            for (auto& elem : this_model_part.Elements()) {
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
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(PLASTIC_STRAIN, aux) - 1.0), tolerance);
                }
            }
        }

        /**
        * Checks the correct work of the internal variable interpolation process after remesh LST test
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(TestInternalInterpolationProcessLST2, KratosInternalInterpolationProcessFastSuite)
        {
            ModelPart this_model_part("Main");
            this_model_part.SetBufferSize(2);
            ProcessInfo& current_process_info = this_model_part.GetProcessInfo();
            current_process_info[DOMAIN_SIZE] = 3;

            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);

            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            auto p_this_law = Kratos::make_shared<LinearJ2Plasticity3D>();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);

            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // First we create the nodes 
            NodeType::Pointer p_node_1 = this_model_part.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_2 = this_model_part.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_3 = this_model_part.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_4 = this_model_part.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_5 = this_model_part.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_6 = this_model_part.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
            
            NodeType::Pointer p_node_7 = this_model_part.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_8 = this_model_part.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
            NodeType::Pointer p_node_9 = this_model_part.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
            NodeType::Pointer p_node_10 = this_model_part.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
            NodeType::Pointer p_node_11 = this_model_part.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
            NodeType::Pointer p_node_12 = this_model_part.CreateNewNode(12 , 2.0 , 0.0 , 0.0);
            
            // Now we create the "conditions"
            std::vector<NodeType::Pointer> element_nodes_0 (4);
            element_nodes_0[0] = p_node_12;
            element_nodes_0[1] = p_node_10;
            element_nodes_0[2] = p_node_8;
            element_nodes_0[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_0( element_nodes_0 );
            
            std::vector<NodeType::Pointer> element_nodes_1 (4);
            element_nodes_1[0] = p_node_4;
            element_nodes_1[1] = p_node_6;
            element_nodes_1[2] = p_node_9;
            element_nodes_1[3] = p_node_7;
            Tetrahedra3D4 <NodeType> tetrahedra_1( element_nodes_1 );
            
            std::vector<NodeType::Pointer> element_nodes_2 (4);
            element_nodes_2[0] = p_node_11;
            element_nodes_2[1] = p_node_7;
            element_nodes_2[2] = p_node_9;
            element_nodes_2[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_2( element_nodes_2 );
            
            std::vector<NodeType::Pointer> element_nodes_3 (4);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_3;
            element_nodes_3[2] = p_node_8;
            element_nodes_3[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_3( element_nodes_3 );
            
            std::vector<NodeType::Pointer> element_nodes_4 (4);
            element_nodes_4[0] = p_node_4;
            element_nodes_4[1] = p_node_6;
            element_nodes_4[2] = p_node_7;
            element_nodes_4[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_4( element_nodes_4 );
            
            std::vector<NodeType::Pointer> element_nodes_5 (4);
            element_nodes_5[0] = p_node_2;
            element_nodes_5[1] = p_node_3;
            element_nodes_5[2] = p_node_5;
            element_nodes_5[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_5( element_nodes_5 );
            
            std::vector<NodeType::Pointer> element_nodes_6 (4);
            element_nodes_6[0] = p_node_10;
            element_nodes_6[1] = p_node_9;
            element_nodes_6[2] = p_node_6;
            element_nodes_6[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_6( element_nodes_6 );
            
            std::vector<NodeType::Pointer> element_nodes_7 (4);
            element_nodes_7[0] = p_node_7;
            element_nodes_7[1] = p_node_8;
            element_nodes_7[2] = p_node_3;
            element_nodes_7[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_7( element_nodes_7 );
            
            std::vector<NodeType::Pointer> element_nodes_8 (4);
            element_nodes_8[0] = p_node_7;
            element_nodes_8[1] = p_node_8;
            element_nodes_8[2] = p_node_6;
            element_nodes_8[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_8( element_nodes_8 );
            
            std::vector<NodeType::Pointer> element_nodes_9 (4);
            element_nodes_9[0] = p_node_4;
            element_nodes_9[1] = p_node_1;
            element_nodes_9[2] = p_node_6;
            element_nodes_9[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_9( element_nodes_9 );
            
            std::vector<NodeType::Pointer> element_nodes_10 (4);
            element_nodes_10[0] = p_node_9;
            element_nodes_10[1] = p_node_12;
            element_nodes_10[2] = p_node_11;
            element_nodes_10[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_10( element_nodes_10 );
            
            std::vector<NodeType::Pointer> element_nodes_11 (4);
            element_nodes_11[0] = p_node_3;
            element_nodes_11[1] = p_node_2;
            element_nodes_11[2] = p_node_1;
            element_nodes_11[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_11( element_nodes_11 );

            // Creating pointers
            auto p_tetrahedra_0 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_0);
            auto p_tetrahedra_1 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_1);
            auto p_tetrahedra_2 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_2);
            auto p_tetrahedra_3 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_3);
            auto p_tetrahedra_4 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_4);
            auto p_tetrahedra_5 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_5);
            auto p_tetrahedra_6 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_6);
            auto p_tetrahedra_7 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_7);
            auto p_tetrahedra_8 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_8);
            auto p_tetrahedra_9 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_9);
            auto p_tetrahedra_10 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_10);
            auto p_tetrahedra_11 = Kratos::make_shared<Tetrahedra3D4 <NodeType>> (tetrahedra_11);

            Element::Pointer p_elem_0 = Kratos::make_shared<SmallDisplacementBbar>(1, p_tetrahedra_0, p_elem_prop);
            Element::Pointer p_elem_1 = Kratos::make_shared<SmallDisplacementBbar>(2, p_tetrahedra_1, p_elem_prop);
            Element::Pointer p_elem_2 = Kratos::make_shared<SmallDisplacementBbar>(3, p_tetrahedra_2, p_elem_prop);
            Element::Pointer p_elem_3 = Kratos::make_shared<SmallDisplacementBbar>(4, p_tetrahedra_3, p_elem_prop);
            Element::Pointer p_elem_4 = Kratos::make_shared<SmallDisplacementBbar>(5, p_tetrahedra_4, p_elem_prop);
            Element::Pointer p_elem_5 = Kratos::make_shared<SmallDisplacementBbar>(6, p_tetrahedra_5, p_elem_prop);
            Element::Pointer p_elem_6 = Kratos::make_shared<SmallDisplacementBbar>(7, p_tetrahedra_6, p_elem_prop);
            Element::Pointer p_elem_7 = Kratos::make_shared<SmallDisplacementBbar>(8, p_tetrahedra_7, p_elem_prop);
            Element::Pointer p_elem_8 = Kratos::make_shared<SmallDisplacementBbar>(9, p_tetrahedra_8, p_elem_prop);
            Element::Pointer p_elem_9 = Kratos::make_shared<SmallDisplacementBbar>(10, p_tetrahedra_9, p_elem_prop);
            Element::Pointer p_elem_10 = Kratos::make_shared<SmallDisplacementBbar>(11, p_tetrahedra_10, p_elem_prop);
            Element::Pointer p_elem_11 = Kratos::make_shared<SmallDisplacementBbar>(12, p_tetrahedra_11, p_elem_prop);
            
            // Adding to the model part
            this_model_part.AddElement(p_elem_0);
            this_model_part.AddElement(p_elem_1);
            this_model_part.AddElement(p_elem_2);
            this_model_part.AddElement(p_elem_3);
            this_model_part.AddElement(p_elem_4);
            this_model_part.AddElement(p_elem_5);
            this_model_part.AddElement(p_elem_6);
            this_model_part.AddElement(p_elem_7);
            this_model_part.AddElement(p_elem_8);
            this_model_part.AddElement(p_elem_9);
            this_model_part.AddElement(p_elem_10);
            this_model_part.AddElement(p_elem_11);

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
            
            // Set DISTANCE and other variables
            Vector ref_metric = ZeroVector(6);
            ref_metric[0] = 1.0;
            ref_metric[3] = 1.0;
            ref_metric[5] = 1.0;
            for (std::size_t i_node = 0; i_node < this_model_part.Nodes().size(); ++i_node) {
                auto it_node = this_model_part.Nodes().begin() + i_node;
                it_node->SetValue(MMG_METRIC, ref_metric);
            }

            // Set PLASTIC_STRAIN on the GP
            for (auto& elem : this_model_part.Elements()) {
                auto& r_this_geometry = elem.GetGeometry();

                // Getting the integration points
                auto this_integration_method = elem.GetIntegrationMethod();
                const auto& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                elem.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

                for (std::size_t i = 0; i <integration_points_number; i++) {
                    constitutive_law_vector[i]->SetValue(PLASTIC_STRAIN, 1.0, current_process_info);
                }
            }

            // Compute remesh
            Parameters params = Parameters(R"({
                "framework"                            : "Lagrangian",
                "internal_variables_parameters"        :
                {
                    "interpolation_type"                   : "LST",
                    "internal_variable_interpolation_list" : ["PLASTIC_STRAIN"]
                },
                "echo_level" : 0
                })" );

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "pre2");

            // Compute NodalH
            FindNodalHProcess process = FindNodalHProcess(this_model_part);
            process.Execute();

            MmgProcess<3> mmg_process = MmgProcess<3>(this_model_part, params);
            mmg_process.Execute();

//             // DEBUG
//             GiDIODebugInternalInterpolation(this_model_part, "2");

            const double tolerance = 1.0e-4;
            for (auto& elem : this_model_part.Elements()) {
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
                    KRATOS_CHECK_LESS_EQUAL(std::abs(constitutive_law_vector[i]->GetValue(PLASTIC_STRAIN, aux) - 1.0), tolerance);
                }
            }
        }
    } // namespace Testing
}  // namespace Kratos.
#endif
