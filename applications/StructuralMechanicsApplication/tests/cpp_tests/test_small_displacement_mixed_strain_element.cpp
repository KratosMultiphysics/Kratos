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
        p_node_1->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_2->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.01;
        p_node_3->FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = 0.02;


        // Compute RHS and LHS
        Vector RHS = ZeroVector(9);
        Matrix LHS = ZeroMatrix(9,9);

        p_element->Initialize(); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_RHS({-49572.6,-49572.6,0.00640046,38461.5,11111.1,0.00177083,11111.1,38461.5,-0.00250463});
        const std::vector<double> expected_LHS_row_0({-1.06838e+06,-299145,918803,384615,-85470.1,918803,683761,384615,918803});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

    void LinearPerturbationField(
        ModelPart &rModelPart,
        const double c_1,
        const double c_2)
    {
        array_1d<double, 3> aux_disp;
        for (auto &r_node : rModelPart.Nodes()) {
            aux_disp[0] = c_1 * r_node.X();
            aux_disp[1] = c_2 * r_node.Y();
            aux_disp[2] = 0.0;
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
            r_node.FastGetSolutionStepValue(VOLUMETRIC_STRAIN) = c_1 + c_2;
        }
    }

    /**
    * Checks the Small Displacement Mixed Strain Element
    * Simple test with known analytical solution
    */
    KRATOS_TEST_CASE_IN_SUITE(SmallDisplacementMixedStrainElementVolumetricDeformationField, KratosStructuralMechanicsFastSuite)
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

        // Initialize the element to initialize the constitutive law
        p_element->Initialize();

        // Set a fake displacement and volumetric strain field to compute the residual
        const double alpha = 1.0;
        const double beta = 1.0;
        LinearPerturbationField(r_model_part, alpha, beta);

        Vector RHS = ZeroVector(9);
        Matrix LHS = ZeroMatrix(9,9);
        p_element->CalculateLocalSystem(LHS, RHS, r_model_part.GetProcessInfo());

        // Perturb the previous displacement and volumetric strain field to compute the residual
        const double alpha_perturbed = 1.25;
        const double beta_perturbed = 1.25;
        LinearPerturbationField(r_model_part, alpha_perturbed, beta_perturbed);

        Vector RHS_perturbed = ZeroVector(9);
        p_element->CalculateRightHandSide(RHS_perturbed, r_model_part.GetProcessInfo());

        // Calculate the perturbation RHS
        const double delta_alpha = alpha_perturbed - alpha;
        const double delta_beta = beta_perturbed - beta;
        LinearPerturbationField(r_model_part, delta_alpha, delta_beta);

        Vector RHS_delta = ZeroVector(9);
        p_element->CalculateRightHandSide(RHS_delta, r_model_part.GetProcessInfo());

        // Check the error
        const Vector RHS_error = RHS_perturbed - (RHS + RHS_delta);

        // Check the LHS
        array_1d<double, 9> perturbation_vector;
        for (auto &r_node : r_model_part.Nodes()) {
            perturbation_vector[(r_node.Id() - 1) * 3] = delta_alpha * r_node.X();
            perturbation_vector[(r_node.Id() - 1) * 3 + 1] = delta_beta * r_node.Y();
            perturbation_vector[(r_node.Id() - 1) * 3 + 2] = delta_alpha + delta_beta;
        }
        const Vector RHS_from_LHS = RHS  - prod(LHS,perturbation_vector);
        const Vector RHS_from_LHS_error = RHS_perturbed - RHS_from_LHS;



        // const Vector RHS_error = RHS_perturbed - (RHS + RHS_perturbation);

        // // Perturb the previous displacement and volumetric strain field to compute the residual
        // array_1d<double, 9> aux_perturbation_vol;

        // Vector RHS_perturbation_vol = -prod(LHS, aux_perturbation_vol);

        // const Vector RHS_error_vol = RHS_perturbed - (RHS + RHS_perturbation_vol);

        // // Perturb the previous displacement and volumetric strain field to compute the residual
        // array_1d<double, 9> aux_perturbation_dev;
        // for (auto &r_node : r_model_part.Nodes()) {
        //     aux_perturbation_dev[(r_node.Id() - 1) * 3] = (alpha_perturbed - alpha) * r_node.X();
        //     aux_perturbation_dev[(r_node.Id() - 1) * 3 + 1] = (beta_perturbed - beta) * r_node.Y();
        //     aux_perturbation_dev[(r_node.Id() - 1) * 3 + 2] = 0.0;
        // }

        // Vector RHS_perturbation_dev = -prod(LHS, aux_perturbation_dev);

        // const Vector RHS_error_dev = RHS_perturbed - (RHS + RHS_perturbation_dev);

        KRATOS_WATCH("")
        KRATOS_WATCH(RHS)
        KRATOS_WATCH(RHS_perturbed)
        KRATOS_WATCH(RHS_delta)
        KRATOS_WATCH(RHS_error)
        KRATOS_WATCH("")
        KRATOS_WATCH(RHS_from_LHS)
        KRATOS_WATCH(RHS_from_LHS_error)
        // KRATOS_WATCH(RHS_perturbation_vol)
        // KRATOS_WATCH(RHS_perturbation_dev)
        // KRATOS_WATCH("")
        // KRATOS_WATCH(RHS_error)
        // KRATOS_WATCH(RHS_error_vol)
        // KRATOS_WATCH(RHS_error_dev)


        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        // const std::vector<double> expected_RHS({-49572.6,-49572.6,0.00640046,38461.5,11111.1,0.00177083,11111.1,38461.5,-0.00250463});
        // const std::vector<double> expected_LHS_row_0({-1.06838e+06,-299145,918803,384615,-85470.1,918803,683761,384615,918803});
        // KRATOS_CHECK_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
        // KRATOS_CHECK_VECTOR_RELATIVE_NEAR(row(LHS,0), expected_LHS_row_0, tolerance)
    }

} // namespace Testing
} // namespace Kratos.
