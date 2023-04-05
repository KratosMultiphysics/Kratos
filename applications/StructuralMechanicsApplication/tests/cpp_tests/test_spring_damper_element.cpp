// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "structural_mechanics_application_variables.h"

#include "custom_elements/spring_damper_element.hpp"

namespace Kratos
{
namespace Testing
{

    void AddDofsElement2D(ModelPart& rModelPart){
        for (auto& r_node : rModelPart.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(ROTATION_Z);
        }
    }

    void AddDofsElement3D(ModelPart& rModelPart) {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);

            r_node.AddDof(ROTATION_X);
            r_node.AddDof(ROTATION_Y);
            r_node.AddDof(ROTATION_Z);
        }
    }


    intrusive_ptr<Element> SetUpElement2D(ModelPart& rModelPart)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(ROTATION);


        // Set the element properties
        const auto p_elem_prop = r_model_part.CreateNewProperties(0);

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 0.1, 0.0, 0.0);

        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{ 5.0, 6.0, 0 };
        p_node_1->FastGetSolutionStepValue(ROTATION) = array_1d<double, 3>{ 0, 0, 7.0 };

        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{ 0.0, 0.0, 0 };
        p_node_2->FastGetSolutionStepValue(ROTATION) = array_1d<double, 3>{ 0, 0, 0.0 };

        AddDofsElement2D(r_model_part);

        const std::vector<ModelPart::IndexType> element_nodes{ 1,2 };
        const auto p_element = r_model_part.CreateNewElement("SpringDamperElement2D2N", 1, element_nodes, p_elem_prop);

        p_element->SetValue(NODAL_DISPLACEMENT_STIFFNESS, array_1d<double, 3>{ 2.0, 3.0, 0 });
        p_element->SetValue(NODAL_ROTATIONAL_STIFFNESS, array_1d<double, 3>{ 0.0, 0.0, 4.0 });

        p_element->SetValue(NODAL_DAMPING_RATIO, array_1d<double, 3>{ 8.0, 9.0, 0 });
        p_element->SetValue(NODAL_ROTATIONAL_DAMPING_RATIO, array_1d<double, 3>{ 0.0, 0.0, 10.0 });

        return p_element;
    }

    // Tests the lhs matrix and the rhs vector of the SpringDamperElement2D2N
    KRATOS_TEST_CASE_IN_SUITE(SpringDamperLocalSystem2D2N, KratosStructuralMechanicsFastSuite)
    {
        // create model part
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        // set up 2d element
        auto p_element = SetUpElement2D(r_model_part);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        // initialize element
        p_element->Initialize(r_process_info); 

        // check element
        p_element->Check(r_process_info);


        // Run test
        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int number_of_dofs = number_of_nodes * 3;
        Matrix calculated_lhs = ZeroMatrix(number_of_dofs,number_of_dofs);
        Vector calculated_rhs = ZeroVector(number_of_dofs);
        p_element->CalculateLocalSystem(calculated_lhs, calculated_rhs, r_process_info);

        // Set expected damping matrix
        Matrix expected_lhs = ZeroMatrix(number_of_dofs);
        expected_lhs(0, 0) = 2.0;
        expected_lhs(0, 3) = -2.0;
        expected_lhs(1, 1) = 3.0;
        expected_lhs(1, 4) = -3.0;
        expected_lhs(2, 2) = 4.0;
        expected_lhs(2, 5) = -4.0;

        expected_lhs(3, 0) = -2.0;
        expected_lhs(3, 3) = 2.0;
        expected_lhs(4, 1) = -3.0;
        expected_lhs(4, 4) = 3.0;
        expected_lhs(5, 2) = -4.0;
        expected_lhs(5, 5) = 4.0;

        // Assert Lhs
        KRATOS_CHECK_MATRIX_NEAR(expected_lhs, calculated_lhs, 1e-10);

        // Set expected rhs
        Vector expected_rhs = ZeroVector(number_of_dofs);
        expected_rhs(0) = - 2.0 * 5.0;
        expected_rhs(1) = - 3.0 * 6.0;
        expected_rhs(2) = - 4.0 * 7.0;
        expected_rhs(3) = 2.0 * 5.0;
        expected_rhs(4) = 3.0 * 6.0;
        expected_rhs(5) = 4.0 * 7.0;

        // Assert rhs
        KRATOS_CHECK_VECTOR_NEAR(expected_rhs, calculated_rhs, 1e-10);
    }


    // Tests the lhs matrix and the rhs vector of the SpringDamperElement2D2N
    KRATOS_TEST_CASE_IN_SUITE(SpringDamperDampingMatrix2D2N, KratosStructuralMechanicsFastSuite)
    {
        // create model part
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

        // set up 2d element
        auto p_element = SetUpElement2D(r_model_part);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        // initialize element
        p_element->Initialize(r_process_info);

        // check element
        p_element->Check(r_process_info);

        // Run test
        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int number_of_dofs = number_of_nodes * 3;
        Matrix calculated_damping_matrix = ZeroMatrix(number_of_dofs, number_of_dofs);
        p_element->CalculateDampingMatrix(calculated_damping_matrix,  r_process_info);

        // Set expected damping matrix
        Matrix expected_damping_matrix = ZeroMatrix(number_of_dofs);
        expected_damping_matrix(0, 0) = 8.0;
        expected_damping_matrix(0, 3) = -8.0;
        expected_damping_matrix(1, 1) = 9.0;
        expected_damping_matrix(1, 4) = -9.0;
        expected_damping_matrix(2, 2) = 10.0;
        expected_damping_matrix(2, 5) = -10.0;

        expected_damping_matrix(3, 0) = -8.0;
        expected_damping_matrix(3, 3) = 8.0;
        expected_damping_matrix(4, 1) = -9.0;
        expected_damping_matrix(4, 4) = 9.0;
        expected_damping_matrix(5, 2) = -10.0;
        expected_damping_matrix(5, 5) = 10.0;

        // Assert
        KRATOS_CHECK_MATRIX_NEAR(expected_damping_matrix, calculated_damping_matrix, 1e-10);

    }
}
}
