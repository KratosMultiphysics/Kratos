// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/gid_io.h"
#include "custom_elements/small_displacement_mixed_strain_element.h"

namespace Kratos
{
namespace Testing
{
    /**
    * Checks the Small Displacement Mixed Strain Element
    * Simple test
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedStrainElement, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUMETRIC_STRAIN);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 0.0 , 1.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 0.0 , 0.0);
        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("SmallDisplacementMixedStrainElement2D3N", 1, element_nodes, p_elem_prop);

        // Set a fake displacement and volumetric strain field to compute the residual
        array_1d<double, 3> aux_disp;
        aux_disp[0] = 0.0;
        aux_disp[1] = 0.0;
        aux_disp[2] = 0.0;
        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
        aux_disp[1] = 0.1;
        p_node_3->FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;

        // Compute RHS and LHS
        Vector RHS = ZeroVector(9);
        Matrix LHS = ZeroMatrix(9,9);

        p_element->Initialize(); // Initialize the element to initialize the constitutive law
        // p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());
        KRATOS_WATCH(RHS)
        p_element->CalculateRightHandSide(RHS, r_model_part.GetProcessInfo());

        // KRATOS_WATCH(LHS)
        KRATOS_WATCH(RHS)

        // KRATOS_CHECK_LESS_EQUAL((0.0223607 - process_info[ERROR_OVERALL])/0.0223607, 1.0e-5);
        // KRATOS_CHECK_LESS_EQUAL((0.148492 - process_info[ENERGY_NORM_OVERALL])/0.148492, 1.0e-5);
    }

} // namespace Testing
} // namespace Kratos.
