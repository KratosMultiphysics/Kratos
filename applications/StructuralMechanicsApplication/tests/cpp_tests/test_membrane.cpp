// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Inigo López
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

#include "custom_elements/membrane_element.hpp"

namespace Kratos
{
namespace Testing
{
    void PrintMatrix(Matrix& rMatrix)
    {
        std::cout.precision(5);
        std::cout << std::scientific;
        std::cout << std::showpos;
        std::cout << std::endl;
        std::cout << std::endl;
        for(unsigned int row = 0; row < rMatrix.size1(); ++row){
            for(unsigned int column = 0; column < rMatrix.size2(); column++){
                if(column == 2 || column == 5){
                    std::cout << " " << rMatrix(row, column) << " |";
                }
                else{
                    std::cout << " " << rMatrix(row, column) << " ";
                }
            }

            std::cout << std::endl;

            if(row ==2|| row == 5){
                for(unsigned int j = 0; j < 14*rMatrix.size1(); j++){
                std::cout << "_" ;
                }
                std::cout << " " << std::endl;
            }
            else{
                for(unsigned int i = 0; i < 3; i++){
                    for(unsigned int j = 0; j < 14*4; j++){
                        std::cout << " " ;
                    }
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }

    void AssignRandomDisplacement(Element::Pointer pElement)
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

    KRATOS_TEST_CASE_IN_SUITE(MembraneElement3D3N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

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
        std::vector<ModelPart::IndexType> element_nodes {1,2,3};
        auto p_element = r_model_part.CreateNewElement("MembraneElement3D3N", 1, element_nodes, p_elem_prop);
        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int dimension = p_element->GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_dofs = number_of_nodes * dimension;

        // Set a fake displacement field to compute the residual
        AssignRandomDisplacement(p_element);
        Vector current_displacement = ZeroVector(dimension*number_of_nodes);
        p_element->GetValuesVector(current_displacement);
        KRATOS_WATCH(current_displacement)

        // Compute RHS and LHS
        Vector RHS_original = ZeroVector(number_of_dofs);
        Matrix LHS_original = ZeroMatrix(number_of_dofs,number_of_dofs);

        p_element->Initialize(); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS_original, RHS_original, r_model_part.GetProcessInfo());

        const double delta = 1e-4;

        // Pinging
        array_1d<double, 3>& disp = p_element->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
        disp[0] += delta;
        p_element->GetValuesVector(current_displacement);
        KRATOS_WATCH(current_displacement)

        Vector RHS_pinged = ZeroVector(number_of_dofs);
        Matrix LHS_pinged = ZeroMatrix(number_of_dofs, number_of_dofs);
        // Compute pinged LHS and RHS
        p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_model_part.GetProcessInfo());

        double sensitivity_fd = -(RHS_pinged(0)-RHS_original(0)) / delta;
        double sensitivity_an = 0.5 * (LHS_original(0,0) + LHS_pinged(0,0));
        double relative_error =  std::abs(sensitivity_fd - sensitivity_an)/std::abs(sensitivity_an)*100.0;
        KRATOS_WATCH(RHS_original)
        KRATOS_WATCH(RHS_pinged)
        KRATOS_WATCH(sensitivity_fd)
        KRATOS_WATCH(sensitivity_an)
        KRATOS_WATCH(relative_error)

        // Unpinging
        disp[0] -= delta;
        p_element->GetValuesVector(current_displacement);
        KRATOS_WATCH(current_displacement)

        PrintMatrix(LHS_original);
        PrintMatrix(LHS_pinged);

        // Vector relative_error = ZeroVector(number_of_dofs*number_of_dofs);
        // for(unsigned int i = 0; i < number_of_nodes; i++){
        //     for(unsigned int j = 0; j < dimension; j++){
        //         // Pinging
        //         array_1d<double, 3>& disp = p_element->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        //         disp[j] += delta;

        //         Vector RHS_pinged = ZeroVector(number_of_dofs);
        //         Matrix LHS_pinged = ZeroMatrix(number_of_dofs, number_of_dofs);
        //         // Compute pinged LHS and RHS
        //         p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_model_part.GetProcessInfo());

        //         // Compute the relative error
        //         for(unsigned int k = 0; k < number_of_dofs; k++){
        //             double sensitivity_fd = -(RHS_pinged(k)-RHS_original(k)) / delta;
        //             double sensitivity_an = 0.5 * (LHS_original(i*dimension + j,k) + LHS_pinged(i*dimension + j,k));
        //             relative_error((i*dimension + j)*number_of_dofs + k) =  std::abs(sensitivity_fd - sensitivity_an)/std::abs(sensitivity_an)*100.0;
        //         }

        //         // Unpinging
        //         disp[j] -= delta;
        //     }
        // }

        // // Check RHS and LHS results
        // const double tolerance = 1.0e-5;
        // const std::vector<double> expected_relative_error({1.46272e-07,1.01309e-10,1.7542e-10,1.38848e-07,5.24272e-11,1.11143e-11,1.546e-07,1.38483e-10,4.36637e-11,9.98607e-11,1.37799e-07,7.96259e-11,4.13558e-11,1.93889e-07,7.57346e-11,5.57454e-11,1.06912e-07,1.22961e-10,1.48198e-10,4.76131e-10,2.56644e-07,7.58107e-11,2.45766e-10,3.09843e-07,1.86351e-11,3.24458e-10,2.19107e-07,6.94375e-08,4.53837e-10,3.21559e-10,9.63228e-08,1.36216e-10,1.80921e-11,1.51967e-10,5.95519e-10,5.44484e-10,3.3914e-10,9.69199e-08,2.88361e-10,4.4912e-11,1.41818e-07,1.03024e-10,2.72377e-10,1.24621e-10,4.7442e-10,6.36844e-10,1.31188e-09,1.54982e-07,3.34607e-10,4.79801e-10,2.51739e-07,6.5135e-10,1.80061e-09,1.33376e-10,7.7201e-08,2.9421e-10,2.87426e-12,2.76943e-10,3.00404e-10,2.12235e-10,1.12454e-07,1.38858e-10,1.98338e-10,1.34145e-10,5.35485e-08,3.15205e-10,5.64514e-11,1.63557e-10,2.06306e-11,1.0259e-10,6.47683e-08,1.26292e-10,1.54817e-10,3.839e-10,1.09561e-07,1.87879e-10,2.5524e-10,1.27795e-10,1.35031e-10,3.60059e-10,1.5045e-07});
        // KRATOS_CHECK_VECTOR_RELATIVE_NEAR(relative_error, expected_relative_error, tolerance)
    }

    KRATOS_TEST_CASE_IN_SUITE(MembraneElement3D4N, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

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
        std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
        auto p_element = r_model_part.CreateNewElement("MembraneElement3D4N", 1, element_nodes, p_elem_prop);

        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int dimension = p_element->GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_dofs = number_of_nodes * dimension;

        // Set a fake displacement field to compute the residual
        //AssignRandomDisplacement(p_element);
        Vector current_displacement = ZeroVector(dimension*number_of_nodes);
        p_element->GetValuesVector(current_displacement);
        KRATOS_WATCH(current_displacement)

        // Compute RHS and LHS
        Vector RHS_original = ZeroVector(number_of_dofs);
        Matrix LHS_original = ZeroMatrix(number_of_dofs,number_of_dofs);

        p_element->Initialize(); // Initialize the element to initialize the constitutive law
        p_element->CalculateLocalSystem(LHS_original, RHS_original, r_model_part.GetProcessInfo());

        const double delta = 1e-4;

        // Pinging
        array_1d<double, 3>& disp = p_element->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
        disp[0] += delta;
        p_element->GetValuesVector(current_displacement);
        KRATOS_WATCH(current_displacement)

        Vector RHS_pinged = ZeroVector(number_of_dofs);
        Matrix LHS_pinged = ZeroMatrix(number_of_dofs, number_of_dofs);
        // Compute pinged LHS and RHS
        p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_model_part.GetProcessInfo());

        double sensitivity_fd = -(RHS_pinged(0)-RHS_original(0)) / delta;
        double sensitivity_an = 0.5 * (LHS_original(0,0) + LHS_pinged(0,0));
        double relative_error =  std::abs(sensitivity_fd - sensitivity_an)/std::abs(sensitivity_an)*100.0;
        KRATOS_WATCH(RHS_original)
        KRATOS_WATCH(RHS_pinged)
        KRATOS_WATCH(sensitivity_fd)
        KRATOS_WATCH(sensitivity_an)
        KRATOS_WATCH(relative_error)

        // Unpinging
        disp[0] -= delta;
        p_element->GetValuesVector(current_displacement);
        KRATOS_WATCH(current_displacement)

        PrintMatrix(LHS_original);
        PrintMatrix(LHS_pinged);

        // Vector relative_error = ZeroVector(number_of_dofs*number_of_dofs);
        // for(unsigned int i = 0; i < number_of_nodes; i++){
        //     for(unsigned int j = 0; j < dimension; j++){
        //         // Pinging
        //         array_1d<double, 3>& disp = p_element->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        //         disp[j] += delta;

        //         Vector RHS_pinged = ZeroVector(number_of_dofs);
        //         Matrix LHS_pinged = ZeroMatrix(number_of_dofs, number_of_dofs);
        //         // Compute pinged LHS and RHS
        //         p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_model_part.GetProcessInfo());

        //         // Compute the relative error
        //         for(unsigned int k = 0; k < number_of_dofs; k++){
        //             double sensitivity_fd = -(RHS_pinged(k)-RHS_original(k)) / delta;
        //             double sensitivity_an = 0.5 * (LHS_original(i*dimension + j,k) + LHS_pinged(i*dimension + j,k));
        //             relative_error((i*dimension + j)*number_of_dofs + k) =  std::abs(sensitivity_fd - sensitivity_an)/std::abs(sensitivity_an)*100.0;
        //         }

        //         // Unpinging
        //         disp[j] -= delta;
        //     }
        // }
        // KRATOS_WATCH(relative_error)

        // Check RHS and LHS results
        // const double tolerance = 1.0e-5;
        // const std::vector<double> expected_relative_error({1.46272e-07,1.01309e-10,1.7542e-10,1.38848e-07,5.24272e-11,1.11143e-11,1.546e-07,1.38483e-10,4.36637e-11,9.98607e-11,1.37799e-07,7.96259e-11,4.13558e-11,1.93889e-07,7.57346e-11,5.57454e-11,1.06912e-07,1.22961e-10,1.48198e-10,4.76131e-10,2.56644e-07,7.58107e-11,2.45766e-10,3.09843e-07,1.86351e-11,3.24458e-10,2.19107e-07,6.94375e-08,4.53837e-10,3.21559e-10,9.63228e-08,1.36216e-10,1.80921e-11,1.51967e-10,5.95519e-10,5.44484e-10,3.3914e-10,9.69199e-08,2.88361e-10,4.4912e-11,1.41818e-07,1.03024e-10,2.72377e-10,1.24621e-10,4.7442e-10,6.36844e-10,1.31188e-09,1.54982e-07,3.34607e-10,4.79801e-10,2.51739e-07,6.5135e-10,1.80061e-09,1.33376e-10,7.7201e-08,2.9421e-10,2.87426e-12,2.76943e-10,3.00404e-10,2.12235e-10,1.12454e-07,1.38858e-10,1.98338e-10,1.34145e-10,5.35485e-08,3.15205e-10,5.64514e-11,1.63557e-10,2.06306e-11,1.0259e-10,6.47683e-08,1.26292e-10,1.54817e-10,3.839e-10,1.09561e-07,1.87879e-10,2.5524e-10,1.27795e-10,1.35031e-10,3.60059e-10,1.5045e-07});
        // KRATOS_CHECK_VECTOR_RELATIVE_NEAR(relative_error, expected_relative_error, tolerance)
    }
}
}