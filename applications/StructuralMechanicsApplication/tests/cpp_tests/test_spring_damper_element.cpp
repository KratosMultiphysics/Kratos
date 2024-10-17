// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

// Project includes
#include "containers/model.h"
#include "structural_mechanics_fast_suite.h"
#include "structural_mechanics_application_variables.h"

#include "custom_elements/spring_damper_element.hpp"

namespace Kratos::Testing {

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

    /**
     * \brief Sets up a 2D spring damper element
     * \param  rModelPart Current model part
     * \return SpringDamper2D element
     */
    auto SetUpElement2D(ModelPart& rModelPart)
    {

        rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(ROTATION);


        // Set the element properties
        const auto p_elem_prop = rModelPart.CreateNewProperties(0);

        // Create the test element
        auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = rModelPart.CreateNewNode(2, 0.1, 0.0, 0.0);

        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{ 5.0, 6.0, 0 };
        p_node_1->FastGetSolutionStepValue(ROTATION) = array_1d<double, 3>{ 0, 0, 7.0 };

        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{ 0.0, 0.0, 0 };
        p_node_2->FastGetSolutionStepValue(ROTATION) = array_1d<double, 3>{ 0, 0, 0.0 };

        AddDofsElement2D(rModelPart);

        const std::vector<ModelPart::IndexType> element_nodes{ 1,2 };
        const auto p_element = rModelPart.CreateNewElement("SpringDamperElement2D", 1, element_nodes, p_elem_prop);

        p_element->SetValue(NODAL_DISPLACEMENT_STIFFNESS, array_1d<double, 3>{ 2.0, 3.0, 0 });
        p_element->SetValue(NODAL_ROTATIONAL_STIFFNESS, array_1d<double, 3>{ 0.0, 0.0, 4.0 });

        p_element->SetValue(NODAL_DAMPING_RATIO, array_1d<double, 3>{ 8.0, 9.0, 0 });
        p_element->SetValue(NODAL_ROTATIONAL_DAMPING_RATIO, array_1d<double, 3>{ 0.0, 0.0, 10.0 });

        return p_element;
    }

    /**
     * \brief Sets up a 3D spring damper element
     * \param  rModelPart Current model part 
     * \return SpringDamper3D element
     */
    auto SetUpElement3D(ModelPart& rModelPart)
    {

        rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(ROTATION);


        // Set the element properties
        const auto p_elem_prop = rModelPart.CreateNewProperties(0);

        // Create the test element
        auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = rModelPart.CreateNewNode(2, 0.1, 0.0, 0.0);

        p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{ 5.0, 6.0, 7.0 };
        p_node_1->FastGetSolutionStepValue(ROTATION) = array_1d<double, 3>{ 8.0, 9.0, 10.0 };

        p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{ 0.0, 0.0, 0 };
        p_node_2->FastGetSolutionStepValue(ROTATION) = array_1d<double, 3>{ 0, 0, 0.0 };

        AddDofsElement3D(rModelPart);

        const std::vector<ModelPart::IndexType> element_nodes{ 1,2 };
        const auto p_element = rModelPart.CreateNewElement("SpringDamperElement3D", 1, element_nodes, p_elem_prop);

        p_element->SetValue(NODAL_DISPLACEMENT_STIFFNESS, array_1d<double, 3>{ 2.0, 3.0, 4.0 });
        p_element->SetValue(NODAL_ROTATIONAL_STIFFNESS, array_1d<double, 3>{ 11.0, 12.0, 13.0 });

        p_element->SetValue(NODAL_DAMPING_RATIO, array_1d<double, 3>{ 14.0, 15.0, 16.0 });
        p_element->SetValue(NODAL_ROTATIONAL_DAMPING_RATIO, array_1d<double, 3>{ 17.0, 18.0, 19.0 });

        return p_element;
    }

    // Tests the lhs matrix and the rhs vector of the SpringDamperElement2D
    KRATOS_TEST_CASE_IN_SUITE(SpringDamperLocalSystem2D, KratosStructuralMechanicsFastSuite)
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
        KRATOS_EXPECT_MATRIX_NEAR(expected_lhs, calculated_lhs, 1e-10);

        // Set expected rhs
        Vector expected_rhs = ZeroVector(number_of_dofs);
        expected_rhs(0) = - 2.0 * 5.0;
        expected_rhs(1) = - 3.0 * 6.0;
        expected_rhs(2) = - 4.0 * 7.0;
        expected_rhs(3) = 2.0 * 5.0;
        expected_rhs(4) = 3.0 * 6.0;
        expected_rhs(5) = 4.0 * 7.0;

        // Assert rhs
        KRATOS_EXPECT_VECTOR_NEAR(expected_rhs, calculated_rhs, 1e-10);
    }


    // Tests the lhs matrix and the rhs vector of the SpringDamperElement2D
    KRATOS_TEST_CASE_IN_SUITE(SpringDamperDampingMatrix2D, KratosStructuralMechanicsFastSuite)
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
        KRATOS_EXPECT_MATRIX_NEAR(expected_damping_matrix, calculated_damping_matrix, 1e-10);

    }


    // Tests the lhs matrix and the rhs vector of the SpringDamperElement3D
    KRATOS_TEST_CASE_IN_SUITE(SpringDamperLocalSystem3D, KratosStructuralMechanicsFastSuite)
    {
        // create model part
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

        // set up 2d element
        auto p_element = SetUpElement3D(r_model_part);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        // initialize element
        p_element->Initialize(r_process_info);

        // check element
        p_element->Check(r_process_info);


        // Run test
        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int number_of_dofs = number_of_nodes * 6;
        Matrix calculated_lhs = ZeroMatrix(number_of_dofs, number_of_dofs);
        Vector calculated_rhs = ZeroVector(number_of_dofs);
        p_element->CalculateLocalSystem(calculated_lhs, calculated_rhs, r_process_info);

        // Set expected damping matrix
        Matrix expected_lhs = ZeroMatrix(number_of_dofs);
        expected_lhs(0, 0) = 2.0;
        expected_lhs(0, 6) = -2.0;
        expected_lhs(1, 1) = 3.0;
        expected_lhs(1, 7) = -3.0;
        expected_lhs(2, 2) = 4.0;
        expected_lhs(2, 8) = -4.0;
        expected_lhs(3, 3) = 11.0;
        expected_lhs(3, 9) = -11.0;
        expected_lhs(4, 4) = 12.0;
        expected_lhs(4, 10) = -12.0;
        expected_lhs(5, 5) = 13.0;
        expected_lhs(5, 11) = -13.0;

        expected_lhs(6, 6) = 2.0;
        expected_lhs(6, 0) = -2.0;
        expected_lhs(7, 7) = 3.0;
        expected_lhs(7, 1) = -3.0;
        expected_lhs(8, 8) = 4.0;
        expected_lhs(8, 2) = -4.0;
        expected_lhs(9, 9) = 11.0;
        expected_lhs(9, 3) = -11.0;
        expected_lhs(10, 10) = 12.0;
        expected_lhs(10, 4) = -12.0;
        expected_lhs(11, 11) = 13.0;
        expected_lhs(11, 5) = -13.0;

        // Assert Lhs
        KRATOS_EXPECT_MATRIX_NEAR(expected_lhs, calculated_lhs, 1e-10);

        // Set expected rhs
        Vector expected_rhs = ZeroVector(number_of_dofs);
        expected_rhs(0) = -2.0 * 5.0;
        expected_rhs(1) = -3.0 * 6.0;
        expected_rhs(2) = -4.0 * 7.0;
        expected_rhs(3) = -11.0 * 8.0;
        expected_rhs(4) = -12.0 * 9.0;
        expected_rhs(5) = -13.0 * 10.0;
        expected_rhs(6) = 2.0 * 5.0;
        expected_rhs(7) = 3.0 * 6.0;
        expected_rhs(8) = 4.0 * 7.0;
        expected_rhs(9) = 11.0 * 8.0;
        expected_rhs(10) = 12.0 * 9.0;
        expected_rhs(11) = 13.0 * 10.0;

        // Assert rhs
        KRATOS_EXPECT_VECTOR_NEAR(expected_rhs, calculated_rhs, 1e-10);
    }



    // Tests the damping matrix of the SpringDamperElement3D
    KRATOS_TEST_CASE_IN_SUITE(SpringDamperDampingMatrix3D, KratosStructuralMechanicsFastSuite)
    {
        // create model part
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);

        // set up 2d element
        auto p_element = SetUpElement3D(r_model_part);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        // initialize element
        p_element->Initialize(r_process_info);

        // check element
        p_element->Check(r_process_info);


        // Run test
        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int number_of_dofs = number_of_nodes * 6;
        Matrix calculated_damping_matrix = ZeroMatrix(number_of_dofs, number_of_dofs);
        p_element->CalculateDampingMatrix(calculated_damping_matrix, r_process_info);

        // Set expected damping matrix
        Matrix expected_damping_matrix = ZeroMatrix(number_of_dofs);
        expected_damping_matrix(0, 0) = 14.0;
        expected_damping_matrix(0, 6) = -14.0;
        expected_damping_matrix(1, 1) = 15.0;
        expected_damping_matrix(1, 7) = -15.0;
        expected_damping_matrix(2, 2) = 16.0;
        expected_damping_matrix(2, 8) = -16.0;
        expected_damping_matrix(3, 3) = 17.0;
        expected_damping_matrix(3, 9) = -17.0;
        expected_damping_matrix(4, 4) = 18.0;
        expected_damping_matrix(4, 10) = -18.0;
        expected_damping_matrix(5, 5) = 19.0;
        expected_damping_matrix(5, 11) = -19.0;

        expected_damping_matrix(6, 6) = 14.0;
        expected_damping_matrix(6, 0) = -14.0;
        expected_damping_matrix(7, 7) = 15.0;
        expected_damping_matrix(7, 1) = -15.0;
        expected_damping_matrix(8, 8) = 16.0;
        expected_damping_matrix(8, 2) = -16.0;
        expected_damping_matrix(9, 9) = 17.0;
        expected_damping_matrix(9, 3) = -17.0;
        expected_damping_matrix(10, 10) = 18.0;
        expected_damping_matrix(10, 4) = -18.0;
        expected_damping_matrix(11, 11) = 19.0;
        expected_damping_matrix(11, 5) = -19.0;

        // Assert damping matrix
        KRATOS_EXPECT_MATRIX_NEAR(expected_damping_matrix, calculated_damping_matrix, 1e-10);

    }
}

