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
            
            Testing::Create2DGeometry(this_model_part, "SmallDisplacementElement2D3N");

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElasticPlaneStrain2DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            p_elem_prop->SetValue(YOUNG_MODULUS, 1.0);
            p_elem_prop->SetValue(POISSON_RATIO, 0.0);

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

            Testing::Create3DGeometry(this_model_part, "SmallDisplacementElement3D4N");

            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            Properties::Pointer p_elem_prop = this_model_part.pGetProperties(0);
            if (!KratosComponents<ConstitutiveLaw>::Has("LinearElastic3DLaw"))
                return void();
            ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElastic3DLaw");
            auto p_this_law = r_clone_cl.Clone();
            p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
            p_elem_prop->SetValue(YOUNG_MODULUS, 1.0);
            p_elem_prop->SetValue(POISSON_RATIO, 0.0);

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
