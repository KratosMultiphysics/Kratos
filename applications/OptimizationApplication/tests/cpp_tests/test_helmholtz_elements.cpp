//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//

// Project includes
#include "containers/model.h"
#include "optimization_application_variables.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(HelmholtzSolidElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_SCALAR);
    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_SCALAR_SOURCE);
    r_model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    // // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);

    // Set flags
    r_model_part.GetProcessInfo().SetValue(COMPUTE_HELMHOLTZ_INVERSE, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_INTEGRATED_FIELD, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_RADIUS, 5.0);

    // // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    for (auto& r_node : r_model_part.Nodes())
        r_node.AddDof(HELMHOLTZ_SCALAR);

    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    auto p_element = r_model_part.CreateNewElement(
        "HelmholtzSolidElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    VariableUtils().SetVariable(HELMHOLTZ_SCALAR, 1.0, r_model_part.Nodes());
    VariableUtils().SetVariable(HELMHOLTZ_SCALAR_SOURCE, 1.0, r_model_part.Nodes());
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes());

    Matrix lhs;
    Vector rhs;
    p_element->CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());

    Matrix ref_lhs(4, 4);
    ref_lhs(0, 0) = 1.251041666667e+01;
    ref_lhs(0, 1) = -4.156250000000e+00;
    ref_lhs(0, 2) = -4.156250000000e+00;
    ref_lhs(0, 3) = -4.156250000000e+00;
    ref_lhs(1, 0) = -4.156250000000e+00;
    ref_lhs(1, 1) = 4.177083333333e+00;
    ref_lhs(1, 2) = 1.041666666667e-02;
    ref_lhs(1, 3) = 1.041666666667e-02;
    ref_lhs(2, 0) = -4.156250000000e+00;
    ref_lhs(2, 1) = 1.041666666667e-02;
    ref_lhs(2, 2) = 4.177083333333e+00;
    ref_lhs(2, 3) = 1.041666666667e-02;
    ref_lhs(3, 0) = -4.156250000000e+00;
    ref_lhs(3, 1) = 1.041666666667e-02;
    ref_lhs(3, 2) = 1.041666666667e-02;
    ref_lhs(3, 3) = 4.177083333333e+00;

    KRATOS_CHECK_MATRIX_NEAR(lhs, ref_lhs, 1e-9);

    Vector ref_rhs(4, -4.166666666667e-02);
    KRATOS_CHECK_VECTOR_NEAR(rhs, ref_rhs, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzVectorSolidElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR);
    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR_SOURCE);
    r_model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    // // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);

    // Set flags
    r_model_part.GetProcessInfo().SetValue(COMPUTE_HELMHOLTZ_INVERSE, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_INTEGRATED_FIELD, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_RADIUS, 5.0);

    // // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(HELMHOLTZ_VECTOR_X);
        r_node.AddDof(HELMHOLTZ_VECTOR_Y);
        r_node.AddDof(HELMHOLTZ_VECTOR_Z);
    }

    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    auto p_element = r_model_part.CreateNewElement(
        "HelmholtzVectorSolidElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    Matrix lhs;
    Vector rhs;
    array_1d<double, 3> one_array;
    one_array[0] = 1.0;
    one_array[1] = 1.0;
    one_array[2] = 1.0;

    VariableUtils().SetVariable(HELMHOLTZ_VECTOR, one_array, r_model_part.Nodes());
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes());
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes());

    p_element->CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());

    Matrix ref_lhs = ZeroMatrix(12, 12);
    ref_lhs(0, 0) = 1.251041666667e+01;
    ref_lhs(0, 3) = -4.156250000000e+00;
    ref_lhs(0, 6) = -4.156250000000e+00;
    ref_lhs(0, 9) = -4.156250000000e+00;
    ref_lhs(1, 1) = 1.251041666667e+01;
    ref_lhs(1, 4) = -4.156250000000e+00;
    ref_lhs(1, 7) = -4.156250000000e+00;
    ref_lhs(1, 10) = -4.156250000000e+00;
    ref_lhs(2, 2) = 1.251041666667e+01;
    ref_lhs(2, 5) = -4.156250000000e+00;
    ref_lhs(2, 8) = -4.156250000000e+00;
    ref_lhs(2, 11) = -4.156250000000e+00;
    ref_lhs(3, 0) = -4.156250000000e+00;
    ref_lhs(3, 3) = 4.177083333333e+00;
    ref_lhs(3, 6) = 1.041666666667e-02;
    ref_lhs(3, 9) = 1.041666666667e-02;
    ref_lhs(4, 1) = -4.156250000000e+00;
    ref_lhs(4, 4) = 4.177083333333e+00;
    ref_lhs(4, 7) = 1.041666666667e-02;
    ref_lhs(4, 10) = 1.041666666667e-02;
    ref_lhs(5, 2) = -4.156250000000e+00;
    ref_lhs(5, 5) = 4.177083333333e+00;
    ref_lhs(5, 8) = 1.041666666667e-02;
    ref_lhs(5, 11) = 1.041666666667e-02;
    ref_lhs(6, 0) = -4.156250000000e+00;
    ref_lhs(6, 3) = 1.041666666667e-02;
    ref_lhs(6, 6) = 4.177083333333e+00;
    ref_lhs(6, 9) = 1.041666666667e-02;
    ref_lhs(7, 1) = -4.156250000000e+00;
    ref_lhs(7, 4) = 1.041666666667e-02;
    ref_lhs(7, 7) = 4.177083333333e+00;
    ref_lhs(7, 10) = 1.041666666667e-02;
    ref_lhs(8, 2) = -4.156250000000e+00;
    ref_lhs(8, 5) = 1.041666666667e-02;
    ref_lhs(8, 8) = 4.177083333333e+00;
    ref_lhs(8, 11) = 1.041666666667e-02;
    ref_lhs(9, 0) = -4.156250000000e+00;
    ref_lhs(9, 3) = 1.041666666667e-02;
    ref_lhs(9, 6) = 1.041666666667e-02;
    ref_lhs(9, 9) = 4.177083333333e+00;
    ref_lhs(10, 1) = -4.156250000000e+00;
    ref_lhs(10, 4) = 1.041666666667e-02;
    ref_lhs(10, 7) = 1.041666666667e-02;
    ref_lhs(10, 10) = 4.177083333333e+00;
    ref_lhs(11, 2) = -4.156250000000e+00;
    ref_lhs(11, 5) = 1.041666666667e-02;
    ref_lhs(11, 8) = 1.041666666667e-02;
    ref_lhs(11, 11) = 4.177083333333e+00;

    KRATOS_CHECK_MATRIX_NEAR(lhs, ref_lhs, 1e-9);

    Vector ref_rhs(12, -4.166666666667e-02);
    KRATOS_CHECK_VECTOR_NEAR(rhs, ref_rhs, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzSurfaceElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_SCALAR);
    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_SCALAR_SOURCE);
    r_model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    // // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);

    // Set flags
    r_model_part.GetProcessInfo().SetValue(COMPUTE_HELMHOLTZ_INVERSE, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_INTEGRATED_FIELD, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_RADIUS, 5.0);

    // // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes())
        r_node.AddDof(HELMHOLTZ_SCALAR);

    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    auto p_element = r_model_part.CreateNewElement(
        "HelmholtzSurfaceElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    VariableUtils().SetVariable(HELMHOLTZ_SCALAR, 1.0, r_model_part.Nodes());
    VariableUtils().SetVariable(HELMHOLTZ_SCALAR_SOURCE, 1.0, r_model_part.Nodes());
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes());

    Matrix lhs;
    Vector rhs;
    p_element->CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());

    Matrix ref_lhs(4, 4);
    ref_lhs(0, 0) = 1.677777777778e+01;
    ref_lhs(0, 1) = -4.111111111111e+00;
    ref_lhs(0, 2) = -8.305555555556e+00;
    ref_lhs(0, 3) = -4.111111111111e+00;
    ref_lhs(1, 0) = -4.111111111111e+00;
    ref_lhs(1, 1) = 1.677777777778e+01;
    ref_lhs(1, 2) = -4.111111111111e+00;
    ref_lhs(1, 3) = -8.305555555556e+00;
    ref_lhs(2, 0) = -8.305555555556e+00;
    ref_lhs(2, 1) = -4.111111111111e+00;
    ref_lhs(2, 2) = 1.677777777778e+01;
    ref_lhs(2, 3) = -4.111111111111e+00;
    ref_lhs(3, 0) = -4.111111111111e+00;
    ref_lhs(3, 1) = -8.305555555556e+00;
    ref_lhs(3, 2) = -4.111111111111e+00;
    ref_lhs(3, 3) = 1.677777777778e+01;

    KRATOS_CHECK_MATRIX_NEAR(lhs, ref_lhs, 1e-9);

    Vector ref_rhs(4, -2.500000000000e-01);
    KRATOS_CHECK_VECTOR_NEAR(rhs, ref_rhs, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzVectorSurfaceElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR);
    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR_SOURCE);
    r_model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    // // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);

    // Set flags
    r_model_part.GetProcessInfo().SetValue(COMPUTE_HELMHOLTZ_INVERSE, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_INTEGRATED_FIELD, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_RADIUS, 5.0);

    // // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(HELMHOLTZ_VECTOR_X);
        r_node.AddDof(HELMHOLTZ_VECTOR_Y);
        r_node.AddDof(HELMHOLTZ_VECTOR_Z);
    }

    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    auto p_element = r_model_part.CreateNewElement(
        "HelmholtzVectorSurfaceElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    array_1d<double, 3> one_array;
    one_array[0] = 1.0;
    one_array[1] = 1.0;
    one_array[2] = 1.0;
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR, one_array, r_model_part.Nodes());
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes());
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes());

    Matrix lhs;
    Vector rhs;
    p_element->CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());

    Matrix ref_lhs = ZeroMatrix(12, 12);
    ref_lhs(0, 0) = 1.677777777778e+01;
    ref_lhs(0, 3) = -4.111111111111e+00;
    ref_lhs(0, 6) = -8.305555555556e+00;
    ref_lhs(0, 9) = -4.111111111111e+00;
    ref_lhs(1, 1) = 1.677777777778e+01;
    ref_lhs(1, 4) = -4.111111111111e+00;
    ref_lhs(1, 7) = -8.305555555556e+00;
    ref_lhs(1, 10) = -4.111111111111e+00;
    ref_lhs(2, 2) = 1.677777777778e+01;
    ref_lhs(2, 5) = -4.111111111111e+00;
    ref_lhs(2, 8) = -8.305555555556e+00;
    ref_lhs(2, 11) = -4.111111111111e+00;
    ref_lhs(3, 0) = -4.111111111111e+00;
    ref_lhs(3, 3) = 1.677777777778e+01;
    ref_lhs(3, 6) = -4.111111111111e+00;
    ref_lhs(3, 9) = -8.305555555556e+00;
    ref_lhs(4, 1) = -4.111111111111e+00;
    ref_lhs(4, 4) = 1.677777777778e+01;
    ref_lhs(4, 7) = -4.111111111111e+00;
    ref_lhs(4, 10) = -8.305555555556e+00;
    ref_lhs(5, 2) = -4.111111111111e+00;
    ref_lhs(5, 5) = 1.677777777778e+01;
    ref_lhs(5, 8) = -4.111111111111e+00;
    ref_lhs(5, 11) = -8.305555555556e+00;
    ref_lhs(6, 0) = -8.305555555556e+00;
    ref_lhs(6, 3) = -4.111111111111e+00;
    ref_lhs(6, 6) = 1.677777777778e+01;
    ref_lhs(6, 9) = -4.111111111111e+00;
    ref_lhs(7, 1) = -8.305555555556e+00;
    ref_lhs(7, 4) = -4.111111111111e+00;
    ref_lhs(7, 7) = 1.677777777778e+01;
    ref_lhs(7, 10) = -4.111111111111e+00;
    ref_lhs(8, 2) = -8.305555555556e+00;
    ref_lhs(8, 5) = -4.111111111111e+00;
    ref_lhs(8, 8) = 1.677777777778e+01;
    ref_lhs(8, 11) = -4.111111111111e+00;
    ref_lhs(9, 0) = -4.111111111111e+00;
    ref_lhs(9, 3) = -8.305555555556e+00;
    ref_lhs(9, 6) = -4.111111111111e+00;
    ref_lhs(9, 9) = 1.677777777778e+01;
    ref_lhs(10, 1) = -4.111111111111e+00;
    ref_lhs(10, 4) = -8.305555555556e+00;
    ref_lhs(10, 7) = -4.111111111111e+00;
    ref_lhs(10, 10) = 1.677777777778e+01;
    ref_lhs(11, 2) = -4.111111111111e+00;
    ref_lhs(11, 5) = -8.305555555556e+00;
    ref_lhs(11, 8) = -4.111111111111e+00;
    ref_lhs(11, 11) = 1.677777777778e+01;

    KRATOS_CHECK_MATRIX_NEAR(lhs, ref_lhs, 1e-9);

    Vector ref_rhs(12, -2.500000000000e-01);
    KRATOS_CHECK_VECTOR_NEAR(rhs, ref_rhs, 1e-9);
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzSolidShapeElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR);
    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR_SOURCE);
    r_model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    // // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    p_elem_prop->SetValue(POISSON_RATIO, 0.3);
    ConstitutiveLaw::Pointer p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get("HelmholtzJacobianStiffened3D").Clone();
    p_elem_prop->SetValue(CONSTITUTIVE_LAW,p_constitutive_law);

    // Set flags
    r_model_part.GetProcessInfo().SetValue(COMPUTE_HELMHOLTZ_INVERSE, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_INTEGRATED_FIELD, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE, 0.16666666);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_RADIUS, 5.0);

    // // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(HELMHOLTZ_VECTOR_X);
        r_node.AddDof(HELMHOLTZ_VECTOR_Y);
        r_node.AddDof(HELMHOLTZ_VECTOR_Z);
    }

    std::vector<ModelPart::IndexType> element_nodes{1, 2, 3, 4};
    auto p_element = r_model_part.CreateNewElement(
        "HelmholtzSolidShapeElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    array_1d<double, 3> one_array;
    one_array[0] = 1.0;
    one_array[1] = 1.0;
    one_array[2] = 1.0;
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR, one_array, r_model_part.Nodes());
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes());
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes());

    Matrix lhs;
    Vector rhs;
    p_element->CalculateLocalSystem(lhs, rhs, r_model_part.GetProcessInfo());

    Matrix ref_lhs = ZeroMatrix(12, 12);
    ref_lhs(0, 0) = 6.917734807692e-02;
    ref_lhs(0, 1) = 2.670940064103e-02;
    ref_lhs(0, 2) = 2.670940064103e-02;
    ref_lhs(0, 3) = -2.697649423077e-02;
    ref_lhs(0, 4) = -1.068376025641e-02;
    ref_lhs(0, 5) = -1.068376025641e-02;
    ref_lhs(0, 6) = -2.670935897436e-04;
    ref_lhs(0, 7) = -1.602564038462e-02;
    ref_lhs(0, 9) = -2.670935897436e-04;
    ref_lhs(0, 11) = -1.602564038462e-02;
    ref_lhs(1, 0) = 2.670940064103e-02;
    ref_lhs(1, 1) = 6.917734807692e-02;
    ref_lhs(1, 2) = 2.670940064103e-02;
    ref_lhs(1, 3) = -1.602564038462e-02;
    ref_lhs(1, 4) = -2.670935897436e-04;
    ref_lhs(1, 6) = -1.068376025641e-02;
    ref_lhs(1, 7) = -2.697649423077e-02;
    ref_lhs(1, 8) = -1.068376025641e-02;
    ref_lhs(1, 10) = -2.670935897436e-04;
    ref_lhs(1, 11) = -1.602564038462e-02;
    ref_lhs(2, 0) = 2.670940064103e-02;
    ref_lhs(2, 1) = 2.670940064103e-02;
    ref_lhs(2, 2) = 6.917734807692e-02;
    ref_lhs(2, 3) = -1.602564038462e-02;
    ref_lhs(2, 5) = -2.670935897436e-04;
    ref_lhs(2, 7) = -1.602564038462e-02;
    ref_lhs(2, 8) = -2.670935897436e-04;
    ref_lhs(2, 9) = -1.068376025641e-02;
    ref_lhs(2, 10) = -1.068376025641e-02;
    ref_lhs(2, 11) = -2.697649423077e-02;
    ref_lhs(3, 0) = -2.697649423077e-02;
    ref_lhs(3, 1) = -1.602564038462e-02;
    ref_lhs(3, 2) = -1.602564038462e-02;
    ref_lhs(3, 3) = 4.780982756410e-02;
    ref_lhs(3, 6) = 1.041666666667e-02;
    ref_lhs(3, 7) = 1.602564038462e-02;
    ref_lhs(3, 9) = 1.041666666667e-02;
    ref_lhs(3, 11) = 1.602564038462e-02;
    ref_lhs(4, 0) = -1.068376025641e-02;
    ref_lhs(4, 1) = -2.670935897436e-04;
    ref_lhs(4, 4) = 2.110042692308e-02;
    ref_lhs(4, 6) = 1.068376025641e-02;
    ref_lhs(4, 7) = 1.041666666667e-02;
    ref_lhs(4, 10) = 1.041666666667e-02;
    ref_lhs(5, 0) = -1.068376025641e-02;
    ref_lhs(5, 2) = -2.670935897436e-04;
    ref_lhs(5, 5) = 2.110042692308e-02;
    ref_lhs(5, 8) = 1.041666666667e-02;
    ref_lhs(5, 9) = 1.068376025641e-02;
    ref_lhs(5, 11) = 1.041666666667e-02;
    ref_lhs(6, 0) = -2.670935897436e-04;
    ref_lhs(6, 1) = -1.068376025641e-02;
    ref_lhs(6, 3) = 1.041666666667e-02;
    ref_lhs(6, 4) = 1.068376025641e-02;
    ref_lhs(6, 6) = 2.110042692308e-02;
    ref_lhs(6, 9) = 1.041666666667e-02;
    ref_lhs(7, 0) = -1.602564038462e-02;
    ref_lhs(7, 1) = -2.697649423077e-02;
    ref_lhs(7, 2) = -1.602564038462e-02;
    ref_lhs(7, 3) = 1.602564038462e-02;
    ref_lhs(7, 4) = 1.041666666667e-02;
    ref_lhs(7, 7) = 4.780982756410e-02;
    ref_lhs(7, 10) = 1.041666666667e-02;
    ref_lhs(7, 11) = 1.602564038462e-02;
    ref_lhs(8, 1) = -1.068376025641e-02;
    ref_lhs(8, 2) = -2.670935897436e-04;
    ref_lhs(8, 5) = 1.041666666667e-02;
    ref_lhs(8, 8) = 2.110042692308e-02;
    ref_lhs(8, 10) = 1.068376025641e-02;
    ref_lhs(8, 11) = 1.041666666667e-02;
    ref_lhs(9, 0) = -2.670935897436e-04;
    ref_lhs(9, 2) = -1.068376025641e-02;
    ref_lhs(9, 3) = 1.041666666667e-02;
    ref_lhs(9, 5) = 1.068376025641e-02;
    ref_lhs(9, 6) = 1.041666666667e-02;
    ref_lhs(9, 9) = 2.110042692308e-02;
    ref_lhs(10, 1) = -2.670935897436e-04;
    ref_lhs(10, 2) = -1.068376025641e-02;
    ref_lhs(10, 4) = 1.041666666667e-02;
    ref_lhs(10, 7) = 1.041666666667e-02;
    ref_lhs(10, 8) = 1.068376025641e-02;
    ref_lhs(10, 10) = 2.110042692308e-02;
    ref_lhs(11, 0) = -1.602564038462e-02;
    ref_lhs(11, 1) = -1.602564038462e-02;
    ref_lhs(11, 2) = -2.697649423077e-02;
    ref_lhs(11, 3) = 1.602564038462e-02;
    ref_lhs(11, 5) = 1.041666666667e-02;
    ref_lhs(11, 7) = 1.602564038462e-02;
    ref_lhs(11, 8) = 1.041666666667e-02;
    ref_lhs(11, 11) = 4.780982756410e-02;

    KRATOS_CHECK_MATRIX_NEAR(lhs, ref_lhs, 1e-9);

    Vector ref_rhs(12, -4.166666666667e-02);
    KRATOS_CHECK_VECTOR_NEAR(rhs, ref_rhs, 1e-9);
}

} // namespace Kratos::Testing
