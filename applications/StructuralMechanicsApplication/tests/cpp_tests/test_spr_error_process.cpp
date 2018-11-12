// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "includes/gid_io.h"

/* Processes */
#include "processes/compute_nodal_gradient_process.h"
#include "custom_processes/spr_error_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        
        void Create2DGeometry(ModelPart& ThisModelPart, const std::string& ElementName)
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
            Triangle2D3 <NodeType> triangle_0( PointerVector<NodeType>{element_nodes_0} );

            std::vector<NodeType::Pointer> element_nodes_1 (3);
            element_nodes_1[0] = p_node_1;
            element_nodes_1[1] = p_node_3;
            element_nodes_1[2] = p_node_4;
            Triangle2D3 <NodeType> triangle_1( PointerVector<NodeType>{element_nodes_1} );

            std::vector<NodeType::Pointer> element_nodes_2 (3);
            element_nodes_2[0] = p_node_2;
            element_nodes_2[1] = p_node_5;
            element_nodes_2[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_2( PointerVector<NodeType>{element_nodes_2} );

            std::vector<NodeType::Pointer> element_nodes_3 (3);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_6;
            element_nodes_3[2] = p_node_3;
            Triangle2D3 <NodeType> triangle_3( PointerVector<NodeType>{element_nodes_3} );

            Element::Pointer p_elem_0 = ThisModelPart.CreateNewElement(ElementName, 1, triangle_0, p_elem_prop);
            Element::Pointer p_elem_1 = ThisModelPart.CreateNewElement(ElementName, 2, triangle_1, p_elem_prop);
            Element::Pointer p_elem_2 = ThisModelPart.CreateNewElement(ElementName, 3, triangle_2, p_elem_prop);
            Element::Pointer p_elem_3 = ThisModelPart.CreateNewElement(ElementName, 4, triangle_3, p_elem_prop);
        }

        void Create3DGeometry(ModelPart& ThisModelPart, const std::string& ElementName)
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
            Tetrahedra3D4 <NodeType> tetrahedra_0( PointerVector<NodeType>{element_nodes_0} );

            std::vector<NodeType::Pointer> element_nodes_1 (4);
            element_nodes_1[0] = p_node_4;
            element_nodes_1[1] = p_node_6;
            element_nodes_1[2] = p_node_9;
            element_nodes_1[3] = p_node_7;
            Tetrahedra3D4 <NodeType> tetrahedra_1( PointerVector<NodeType>{element_nodes_1} );

            std::vector<NodeType::Pointer> element_nodes_2 (4);
            element_nodes_2[0] = p_node_11;
            element_nodes_2[1] = p_node_7;
            element_nodes_2[2] = p_node_9;
            element_nodes_2[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_2( PointerVector<NodeType>{element_nodes_2} );

            std::vector<NodeType::Pointer> element_nodes_3 (4);
            element_nodes_3[0] = p_node_5;
            element_nodes_3[1] = p_node_3;
            element_nodes_3[2] = p_node_8;
            element_nodes_3[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_3( PointerVector<NodeType>{element_nodes_3} );

            std::vector<NodeType::Pointer> element_nodes_4 (4);
            element_nodes_4[0] = p_node_4;
            element_nodes_4[1] = p_node_6;
            element_nodes_4[2] = p_node_7;
            element_nodes_4[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_4( PointerVector<NodeType>{element_nodes_4} );

            std::vector<NodeType::Pointer> element_nodes_5 (4);
            element_nodes_5[0] = p_node_2;
            element_nodes_5[1] = p_node_3;
            element_nodes_5[2] = p_node_5;
            element_nodes_5[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_5( PointerVector<NodeType>{element_nodes_5} );

            std::vector<NodeType::Pointer> element_nodes_6 (4);
            element_nodes_6[0] = p_node_10;
            element_nodes_6[1] = p_node_9;
            element_nodes_6[2] = p_node_6;
            element_nodes_6[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_6( PointerVector<NodeType>{element_nodes_6} );

            std::vector<NodeType::Pointer> element_nodes_7 (4);
            element_nodes_7[0] = p_node_7;
            element_nodes_7[1] = p_node_8;
            element_nodes_7[2] = p_node_3;
            element_nodes_7[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_7( PointerVector<NodeType>{element_nodes_7} );

            std::vector<NodeType::Pointer> element_nodes_8 (4);
            element_nodes_8[0] = p_node_7;
            element_nodes_8[1] = p_node_8;
            element_nodes_8[2] = p_node_6;
            element_nodes_8[3] = p_node_9;
            Tetrahedra3D4 <NodeType> tetrahedra_8( PointerVector<NodeType>{element_nodes_8} );

            std::vector<NodeType::Pointer> element_nodes_9 (4);
            element_nodes_9[0] = p_node_4;
            element_nodes_9[1] = p_node_1;
            element_nodes_9[2] = p_node_6;
            element_nodes_9[3] = p_node_3;
            Tetrahedra3D4 <NodeType> tetrahedra_9( PointerVector<NodeType>{element_nodes_9} );

            std::vector<NodeType::Pointer> element_nodes_10 (4);
            element_nodes_10[0] = p_node_9;
            element_nodes_10[1] = p_node_12;
            element_nodes_10[2] = p_node_11;
            element_nodes_10[3] = p_node_8;
            Tetrahedra3D4 <NodeType> tetrahedra_10( PointerVector<NodeType>{element_nodes_10} );

            std::vector<NodeType::Pointer> element_nodes_11 (4);
            element_nodes_11[0] = p_node_3;
            element_nodes_11[1] = p_node_2;
            element_nodes_11[2] = p_node_1;
            element_nodes_11[3] = p_node_6;
            Tetrahedra3D4 <NodeType> tetrahedra_11( PointerVector<NodeType>{element_nodes_11} );

            Element::Pointer p_elem_0 = ThisModelPart.CreateNewElement(ElementName, 1, tetrahedra_0, p_elem_prop);
            Element::Pointer p_elem_1 = ThisModelPart.CreateNewElement(ElementName, 2, tetrahedra_1, p_elem_prop);
            Element::Pointer p_elem_2 = ThisModelPart.CreateNewElement(ElementName, 3, tetrahedra_2, p_elem_prop);
            Element::Pointer p_elem_3 = ThisModelPart.CreateNewElement(ElementName, 4, tetrahedra_3, p_elem_prop);
            Element::Pointer p_elem_4 = ThisModelPart.CreateNewElement(ElementName, 5, tetrahedra_4, p_elem_prop);
            Element::Pointer p_elem_5 = ThisModelPart.CreateNewElement(ElementName, 6, tetrahedra_5, p_elem_prop);
            Element::Pointer p_elem_6 = ThisModelPart.CreateNewElement(ElementName, 7, tetrahedra_6, p_elem_prop);
            Element::Pointer p_elem_7 = ThisModelPart.CreateNewElement(ElementName, 8, tetrahedra_7, p_elem_prop);
            Element::Pointer p_elem_8 = ThisModelPart.CreateNewElement(ElementName, 9, tetrahedra_8, p_elem_prop);
            Element::Pointer p_elem_9 = ThisModelPart.CreateNewElement(ElementName, 10, tetrahedra_9, p_elem_prop);
            Element::Pointer p_elem_10 = ThisModelPart.CreateNewElement(ElementName, 11, tetrahedra_10, p_elem_prop);
            Element::Pointer p_elem_11 = ThisModelPart.CreateNewElement(ElementName, 12, tetrahedra_11, p_elem_prop);
        }

        /**
        * Checks the correct work of the SPR metric process
        * Test triangle 
        */

        KRATOS_TEST_CASE_IN_SUITE(SPRErrorProcess1, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main",2);
            
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);
            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElasticPlaneStrain2DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            p_elem_prop->SetValue(YOUNG_MODULUS, 1.0);
            p_elem_prop->SetValue(POISSON_RATIO, 0.0);

            Testing::Create2DGeometry(this_model_part, "SmallDisplacementElement2D3N");

            for (auto& node : this_model_part.Nodes()) {
                if (node.X() > 0.9) {
                    node.FastGetSolutionStepValue(DISPLACEMENT_X) += 0.1;
                    node.Coordinates()[0] += 0.1;
                }
            }

            for (auto& ielem : this_model_part.Elements()) {
                ielem.Initialize();
                ielem.InitializeSolutionStep(process_info);
            }

            // Compute error
            SPRErrorProcess<2> spr_process = SPRErrorProcess<2>(this_model_part);
            spr_process.Execute();

            KRATOS_CHECK_LESS_EQUAL((0.0223607 - process_info[ERROR_OVERALL])/0.0223607, 1.0e-5);
            KRATOS_CHECK_LESS_EQUAL((0.148492 - process_info[ENERGY_NORM_OVERALL])/0.148492, 1.0e-5);
        }
        
        /** 
        * Checks the correct work of the nodal SPR compute
        * Test tetrahedra
        */

        KRATOS_TEST_CASE_IN_SUITE(SPRErrorProcess2, KratosStructuralMechanicsFastSuite)
        {
            Model current_model;
            ModelPart& this_model_part = current_model.CreateModelPart("Main",2);
            
            this_model_part.AddNodalSolutionStepVariable(NODAL_H);
            this_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            
            auto& process_info = this_model_part.GetProcessInfo();
            process_info[STEP] = 1;
            process_info[NL_ITERATION_NUMBER] = 1;
            
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            p_elem_prop->SetValue(YOUNG_MODULUS, 1.0);
            p_elem_prop->SetValue(POISSON_RATIO, 0.0);

            Testing::Create3DGeometry(this_model_part, "SmallDisplacementElement3D4N");

            for (auto& node : this_model_part.Nodes()) {
                if (node.X() > 0.9) {
                    node.FastGetSolutionStepValue(DISPLACEMENT_X) += 0.1;
                    node.Coordinates()[0] += 0.1;
                }
            }

            for (auto& ielem : this_model_part.Elements()) {
                ielem.Initialize();
                ielem.InitializeSolutionStep(process_info);
            }

            // Compute error
            SPRErrorProcess<3> spr_process = SPRErrorProcess<3>(this_model_part);
            spr_process.Execute();

            KRATOS_CHECK_LESS_EQUAL((0.122409 - process_info[ERROR_OVERALL])/0.122409, 1.0e-5);
            KRATOS_CHECK_LESS_EQUAL((0.257196 - process_info[ENERGY_NORM_OVERALL])/0.257196, 1.0e-5);
        }
    } // namespace Testing
}  // namespace Kratos.
