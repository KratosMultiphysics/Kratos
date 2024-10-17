// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/line_2d_2.h"
#include "structural_mechanics_fast_suite.h"

// Application includes
#include "custom_elements/membrane_element_2D2N.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MembraneElement2D2N, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("TestModelPart",1);

    // Add required variables
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    Vector prestress(1);
    prestress[0] = 1.0e3;
    p_elem_prop->SetValue(DENSITY, 1.0e3);
    p_elem_prop->SetValue(THICKNESS, 0.1);
    p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e6);
    p_elem_prop->SetValue(PRESTRESS_VECTOR, prestress);

    // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    std::vector<ModelPart::IndexType> element_nodes {1, 2};
    auto p_element = r_model_part.CreateNewElement("MembraneElement2D2N", 1, element_nodes, p_elem_prop);

    // Set a fake displacement and body force fields to compute the residual
    array_1d<double, 3> aux_vect = ZeroVector(3);
    aux_vect[1] = 0.1;
    p_node_1->FastGetSolutionStepValue(DISPLACEMENT) = aux_vect;
    p_node_2->FastGetSolutionStepValue(DISPLACEMENT) = aux_vect;
    aux_vect[1] = -10.0;
    p_node_1->FastGetSolutionStepValue(VOLUME_ACCELERATION) = aux_vect;
    p_node_2->FastGetSolutionStepValue(VOLUME_ACCELERATION) = aux_vect;

    // Compute RHS, LHS and mass matrix
    Vector RHS = ZeroVector(4);
    Matrix LHS = ZeroMatrix(4,4);
    Matrix mass_matrix = ZeroMatrix(4,4);

    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS, RHS, r_process_info);
    p_element->CalculateMassMatrix(mass_matrix, r_process_info);

    // Check RHS, LHS and mass matrix results
    Matrix expected_LHS = ZeroMatrix(4,4);
    expected_LHS(0,0) = 200100.0; expected_LHS(0,2) = -200100.0;
    expected_LHS(1,1) = 100.0; expected_LHS(1,3) = -100.0;
    expected_LHS(2,0) = -200100.0; expected_LHS(2,2) = 200100.0;
    expected_LHS(3,1) = -100.0; expected_LHS(3,3) = 100.0;

    Matrix expected_mass_matrix = ZeroMatrix(4,4);
    expected_mass_matrix(0,0) = 100.0/3.0; expected_mass_matrix(0,2) = 100.0/6.0;
    expected_mass_matrix(1,1) = 100.0/3.0; expected_mass_matrix(1,3) = 100.0/6.0;
    expected_mass_matrix(2,0) = 100.0/6.0; expected_mass_matrix(2,2) = 100.0/3.0;
    expected_mass_matrix(3,1) = 100.0/6.0; expected_mass_matrix(3,3) = 100.0/3.0;

    const std::vector<double> expected_RHS({100, -500, -100, -500});

    const double tolerance = 1.0e-7;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(RHS, expected_RHS, tolerance)
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(LHS, expected_LHS, tolerance)
    KRATOS_EXPECT_MATRIX_RELATIVE_NEAR(mass_matrix, expected_mass_matrix, tolerance)
}

} // namespace Kratos::Testing
