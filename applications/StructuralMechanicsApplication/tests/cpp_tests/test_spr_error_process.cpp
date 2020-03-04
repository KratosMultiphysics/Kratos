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
#include "testing/testing.h"
#include "includes/kratos_flags.h"
#include "utilities/cpp_tests_utilities.h"

/* Processes */
#include "processes/compute_nodal_gradient_process.h"
#include "custom_processes/spr_error_process.h"

namespace Kratos 
{
    namespace Testing 
    {
        typedef Node<3> NodeType;
        
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
            
            CppTestsUtilities::Create2DGeometry(this_model_part, "SmallDisplacementElement2D3N", false);

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

            KRATOS_CHECK_RELATIVE_NEAR(0.0229129, process_info[ERROR_OVERALL], 1.0e-5);
            KRATOS_CHECK_RELATIVE_NEAR(0.148492, process_info[ENERGY_NORM_OVERALL], 1.0e-5);
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

            CppTestsUtilities::Create3DGeometry(this_model_part, "SmallDisplacementElement3D4N", false);

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

            KRATOS_CHECK_RELATIVE_NEAR(0.0494477, process_info[ERROR_OVERALL], 1.0e-5);
            KRATOS_CHECK_RELATIVE_NEAR(0.257196, process_info[ENERGY_NORM_OVERALL], 1.0e-5);
        }
    } // namespace Testing
}  // namespace Kratos.
