// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter and Inigo Lopez
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

#include "custom_elements/membrane_element.hpp"

namespace Kratos
{
namespace Testing
{

    void AddDisplacementDofs(ModelPart& rModelPart){
        for (auto& r_node : rModelPart.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }
    }

    void AssignPredefinedDisplacement(Element::Pointer pElement)
    {
        const unsigned int number_of_nodes = pElement->GetGeometry().size();
        const unsigned int dimension = pElement->GetGeometry().WorkingSpaceDimension();

        double displacement = 0.0;
        for(unsigned int i = 0; i < number_of_nodes; i++){
            array_1d<double, 3>& disp = pElement->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            for(unsigned int j = 0; j < dimension; j++){
                disp[j] = displacement;
                displacement += 0.1;
            }
        }
    }

    void ComputeRelativeError(Vector& rRelativeError, Element::Pointer pElement, ModelPart& rModelPart)
    {
        const unsigned int number_of_nodes = pElement->GetGeometry().size();
        const unsigned int dimension = pElement->GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_dofs = number_of_nodes * dimension;
        // Set a predifined displacement field to compute the residual
        AssignPredefinedDisplacement(pElement);

        // Compute RHS and LHS
        Vector RHS_original = ZeroVector(number_of_dofs);
        Matrix LHS_original = ZeroMatrix(number_of_dofs,number_of_dofs);

        pElement->Initialize(rModelPart.GetProcessInfo()); // Initialize the element to initialize the constitutive law
        pElement->Check(rModelPart.GetProcessInfo());
        pElement->CalculateLocalSystem(LHS_original, RHS_original, rModelPart.GetProcessInfo());

        const double delta = 1e-4;

        for(unsigned int i = 0; i < number_of_nodes; i++){
            for(unsigned int j = 0; j < dimension; j++){
                // Pinging
                array_1d<double, 3>& disp = pElement->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
                disp[j] += delta;

                Vector RHS_pinged = ZeroVector(number_of_dofs);
                Matrix LHS_pinged = ZeroMatrix(number_of_dofs, number_of_dofs);
                // Compute pinged LHS and RHS
                pElement->CalculateLocalSystem(LHS_pinged, RHS_pinged, rModelPart.GetProcessInfo());

                for(unsigned int k = 0; k < number_of_dofs; k++){
                    // Compute the finite difference estimate of the sensitivity
                    double sensitivity_fd = -(RHS_pinged(k)-RHS_original(k)) / delta;
                    // Compute the average of the original and pinged analytic sensitivities
                    double sensitivity_an = 0.5 * (LHS_original(i*dimension + j,k) + LHS_pinged(i*dimension + j,k));
                    // Compute the relative error between sensitivities
                    rRelativeError((i*dimension + j)*number_of_dofs + k) =  std::abs(sensitivity_fd - sensitivity_an)/std::abs(sensitivity_an)*100.0;
                }

                // Unpinging
                disp[j] -= delta;
            }
        }
    }

    // Tests the stiffness matrix of the MembraneElement3D3N comparing analytical and numerical elemental sensitivities
    KRATOS_TEST_CASE_IN_SUITE(MembraneElement3D3N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);

        AddDisplacementDofs(r_model_part);


        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("MembraneElement3D3N", 1, element_nodes, p_elem_prop);
        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int dimension = p_element->GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_dofs = number_of_nodes * dimension;


        Vector relative_error = ZeroVector(number_of_dofs*number_of_dofs);
        ComputeRelativeError(relative_error, p_element, r_model_part);

        // Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_relative_error({1.46272e-07,1.01309e-10,1.7542e-10,1.38848e-07,5.24272e-11,1.11143e-11,1.546e-07,1.38483e-10,4.36637e-11,9.98607e-11,1.37799e-07,7.96259e-11,4.13558e-11,1.93889e-07,7.57346e-11,5.57454e-11,1.06912e-07,1.22961e-10,1.48198e-10,4.76131e-10,2.56644e-07,7.58107e-11,2.45766e-10,3.09843e-07,1.86351e-11,3.24458e-10,2.19107e-07,6.94375e-08,4.53837e-10,3.21559e-10,9.63228e-08,1.36216e-10,1.80921e-11,1.51967e-10,5.95519e-10,5.44484e-10,3.3914e-10,9.69199e-08,2.88361e-10,4.4912e-11,1.41818e-07,1.03024e-10,2.72377e-10,1.24621e-10,4.7442e-10,6.36844e-10,1.31188e-09,1.54982e-07,3.34607e-10,4.79801e-10,2.51739e-07,6.5135e-10,1.80061e-09,1.33376e-10,7.7201e-08,2.9421e-10,2.87426e-12,2.76943e-10,3.00404e-10,2.12235e-10,1.12454e-07,1.38858e-10,1.98338e-10,1.34145e-10,5.35485e-08,3.15205e-10,5.64514e-11,1.63557e-10,2.06306e-11,1.0259e-10,6.47683e-08,1.26292e-10,1.54817e-10,3.839e-10,1.09561e-07,1.87879e-10,2.5524e-10,1.27795e-10,1.35031e-10,3.60059e-10,1.5045e-07});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(relative_error, expected_relative_error, tolerance)
    }

    // Tests the stiffness matrix of the MembraneElement3D4N comparing analytical and numerical elemental sensitivities
    KRATOS_TEST_CASE_IN_SUITE(MembraneElement3D4N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 1.0 , 0.0);

        AddDisplacementDofs(r_model_part);

        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("MembraneElement3D4N", 1, element_nodes, p_elem_prop);

        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int dimension = p_element->GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_dofs = number_of_nodes * dimension;

        Vector relative_error = ZeroVector(number_of_dofs*number_of_dofs);
        ComputeRelativeError(relative_error, p_element, r_model_part);

        //Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_relative_error({7.83999e-08,2.07838e-11,2.85724e-10,1.19805e-07,2.53715e-09,6.58355e-10,6.34008e-08,5.77027e-11,3.9188e-10,7.30737e-08,2.26601e-10,2.88665e-10,1.56401e-10,6.67553e-08,2.00164e-10,3.7924e-09,4.32224e-07,2.78981e-09,3.98132e-10,5.89574e-08,2.76404e-10,1.40101e-10,3.35175e-08,1.16888e-10,5.6946e-10,5.4851e-10,1.3056e-07,2.87022e-10,2.06904e-12,3.2299e-05,4.79543e-10,9.35669e-11,1.21407e-07,1.03257e-09,2.52929e-10,6.84278e-08,1.20185e-07,1.68722e-08,9.15427e-10,1.88326e-07,3.92703e-10,5.80919e-11,3.79004e-07,7.25958e-10,8.25508e-10,2.00979e-07,1.47245e-09,1.6338e-08,1.42025e-09,4.38993e-07,6.80859e-09,2.92253e-11,1.19192e-07,4.58938e-11,1.20997e-10,7.48042e-08,1.48577e-10,3.63387e-09,7.41823e-08,1.16976e-11,1.08821e-09,8.9436e-09,3.25153e-05,3.11555e-10,4.59753e-11,3.59894e-07,2.29186e-09,6.19311e-10,2.72257e-07,1.00945e-08,1.85752e-09,2.25757e-07,6.36167e-08,2.45672e-10,9.80188e-10,3.78806e-07,1.66561e-10,3.77833e-10,1.51191e-07,1.15385e-10,4.436e-10,8.086e-06,4.61353e-11,1.95153e-10,2.45672e-10,5.88961e-08,1.05936e-09,4.80807e-11,7.45737e-08,1.67813e-10,2.08987e-10,1.12515e-07,5.10433e-11,1.02384e-09,1.9622e-07,7.09616e-10,1.02271e-10,1.0528e-09,1.21355e-07,7.72253e-10,3.54967e-10,2.72712e-07,4.1818e-10,1.81686e-10,3.05914e-07,8.06128e-10,1.81512e-09,3.57271e-07,7.28385e-08,7.78619e-10,3.12537e-10,2.00335e-07,3.25997e-09,4.32728e-09,8.16507e-06,1.47387e-09,7.04346e-10,1.60776e-07,1.28813e-09,1.48714e-09,3.26769e-10,3.37628e-08,6.24457e-10,7.40437e-10,7.46591e-08,9.46029e-11,3.09336e-10,1.9826e-07,9.3412e-10,1.80188e-10,7.87286e-08,2.51837e-10,1.04632e-09,4.97367e-10,6.78623e-08,7.14672e-09,3.06683e-10,2.24474e-07,2.04464e-09,3.84428e-09,3.56122e-07,1.10454e-09,3.60626e-10,1.84727e-07});
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(relative_error, expected_relative_error, tolerance)
    }
}
}
