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
#include "testing/testing.h"
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "optimization_application_variables.h"
#include "includes/debug_helpers.h"

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(HelmholtzSolidElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
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
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0,0.0,0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0,1.0,0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0,0.0,1.0);

    for (auto& r_node : r_model_part.Nodes())
        r_node.AddDof(HELMHOLTZ_SCALAR);

    std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
    auto p_element = r_model_part.CreateNewElement("HelmholtzSolidElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    VariableUtils().SetVariable(HELMHOLTZ_SCALAR, 1.0, r_model_part.Nodes());
    VariableUtils().SetVariable(HELMHOLTZ_SCALAR_SOURCE, 1.0, r_model_part.Nodes());
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes());

    Matrix lhs;
    Vector rhs;
    p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

    Matrix ref_lhs(4, 4);
    ref_lhs(0, 0) =  1.251041666667e+01;
    ref_lhs(0, 1) = -4.156250000000e+00;
    ref_lhs(0, 2) = -4.156250000000e+00;
    ref_lhs(0, 3) = -4.156250000000e+00;
    ref_lhs(1, 0) = -4.156250000000e+00;
    ref_lhs(1, 1) =  4.177083333333e+00;
    ref_lhs(1, 2) =  1.041666666667e-02;
    ref_lhs(1, 3) =  1.041666666667e-02;
    ref_lhs(2, 0) = -4.156250000000e+00;
    ref_lhs(2, 1) =  1.041666666667e-02;
    ref_lhs(2, 2) =  4.177083333333e+00;
    ref_lhs(2, 3) =  1.041666666667e-02;
    ref_lhs(3, 0) = -4.156250000000e+00;
    ref_lhs(3, 1) =  1.041666666667e-02;
    ref_lhs(3, 2) =  1.041666666667e-02;
    ref_lhs(3, 3) =  4.177083333333e+00;

    KRATOS_CHECK_MATRIX_NEAR(lhs, ref_lhs, 1e-9);

    Vector ref_rhs(4);
    ref_rhs[0] = -4.166666666667e-02;
    ref_rhs[1] = -4.166666666667e-02;
    ref_rhs[2] = -4.166666666667e-02;
    ref_rhs[3] = -4.166666666667e-02;

    KRATOS_CHECK_VECTOR_NEAR(rhs, ref_rhs, 1e-9);

    Vector x;
    MathUtils<double>::Solve(lhs,x,rhs);

    const std::vector<double> reference_x = {-1.0,-1.0,-1.0,-1.0};
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzVectorSolidElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
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
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0,0.0,0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0,1.0,0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0,0.0,1.0);

    for (auto& r_node : r_model_part.Nodes()){
        r_node.AddDof(HELMHOLTZ_VECTOR_X);
        r_node.AddDof(HELMHOLTZ_VECTOR_Y);
        r_node.AddDof(HELMHOLTZ_VECTOR_Z);
    }

    std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
    auto p_element = r_model_part.CreateNewElement("HelmholtzVectorSolidElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    Matrix lhs;
    Vector rhs;
    array_1d<double, 3> one_array;
    one_array[0] = 1.0;
    one_array[1] = 1.0;
    one_array[2] = 1.0;

    VariableUtils().SetVariable(HELMHOLTZ_VECTOR, one_array, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes(), ACTIVE, false);

    p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

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

    Vector ref_rhs(12);
    ref_rhs[0] = -4.166666666667e-02;
    ref_rhs[1] = -4.166666666667e-02;
    ref_rhs[2] = -4.166666666667e-02;
    ref_rhs[3] = -4.166666666667e-02;
    ref_rhs[4] = -4.166666666667e-02;
    ref_rhs[5] = -4.166666666667e-02;
    ref_rhs[6] = -4.166666666667e-02;
    ref_rhs[7] = -4.166666666667e-02;
    ref_rhs[8] = -4.166666666667e-02;
    ref_rhs[9] = -4.166666666667e-02;
    ref_rhs[10] = -4.166666666667e-02;
    ref_rhs[11] = -4.166666666667e-02;

    KRATOS_CHECK_VECTOR_NEAR(rhs, ref_rhs, 1e-9);

    Vector x;
    MathUtils<double>::Solve(lhs,x,rhs);

    const std::vector<double> reference_x = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};

    KRATOS_CHECK_VECTOR_NEAR(x, reference_x, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzSurfaceElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
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
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0,0.0,0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 1.0,1.0,0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0,1.0,0.0);

    for (auto& r_node : r_model_part.Nodes())
        r_node.AddDof(HELMHOLTZ_SCALAR);

    std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
    auto p_element = r_model_part.CreateNewElement("HelmholtzSurfaceElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    Matrix lhs;
    Vector rhs;
    VariableUtils().SetVariable(HELMHOLTZ_SCALAR, 1.0, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(HELMHOLTZ_SCALAR_SOURCE, 1.0, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes(), ACTIVE, false);

    p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

    Vector x;
    MathUtils<double>::Solve(lhs,x,rhs);

    const std::vector<double> reference_x = {-1.0,-1.0,-1.0,-1.0};

    KRATOS_CHECK_VECTOR_NEAR(x, reference_x, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzVectorSurfaceElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
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
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0,0.0,0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 1.0,1.0,0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0,1.0,0.0);

    for (auto& r_node : r_model_part.Nodes()){
        r_node.AddDof(HELMHOLTZ_VECTOR_X);
        r_node.AddDof(HELMHOLTZ_VECTOR_Y);
        r_node.AddDof(HELMHOLTZ_VECTOR_Z);
    }

    std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
    auto p_element = r_model_part.CreateNewElement("HelmholtzVectorSurfaceElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    Matrix lhs;
    Vector rhs;
    array_1d<double, 3> one_array;
    one_array[0] = 1.0;
    one_array[1] = 1.0;
    one_array[2] = 1.0;
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR, one_array, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes(), ACTIVE, false);

    p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

    Vector x;
    MathUtils<double>::Solve(lhs,x,rhs);

    const std::vector<double> reference_x = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};

    KRATOS_CHECK_VECTOR_NEAR(x, reference_x, 1.0e-6);
}

KRATOS_TEST_CASE_IN_SUITE(HelmholtzSolidShapeElement, KratosOptimizationFastSuite)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR);
    r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR_SOURCE);
    r_model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS);

    // // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);

    // Set flags
    r_model_part.GetProcessInfo().SetValue(COMPUTE_HELMHOLTZ_INVERSE, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_INTEGRATED_FIELD, false);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE, 0.16666666);
    r_model_part.GetProcessInfo().SetValue(HELMHOLTZ_RADIUS, 5.0);

    // // Create the test element
    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0,0.0,0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 1.0,0.0,0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0,1.0,0.0);
    auto p_node_4 = r_model_part.CreateNewNode(4, 0.0,0.0,1.0);

    for (auto& r_node : r_model_part.Nodes()){
        r_node.AddDof(HELMHOLTZ_VECTOR_X);
        r_node.AddDof(HELMHOLTZ_VECTOR_Y);
        r_node.AddDof(HELMHOLTZ_VECTOR_Z);
    }

    std::vector<ModelPart::IndexType> element_nodes {1,2,3,4};
    auto p_element = r_model_part.CreateNewElement("HelmholtzSolidShapeElement3D4N", 1, element_nodes, p_elem_prop);

    p_element->Initialize(r_model_part.GetProcessInfo());

    Matrix lhs;
    Vector rhs;
    array_1d<double, 3> one_array;
    one_array[0] = 1.0;
    one_array[1] = 1.0;
    one_array[2] = 1.0;
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR, one_array, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes(), ACTIVE, false);
    VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes(), ACTIVE, false);

    p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

    Vector x;
    MathUtils<double>::Solve(lhs,x,rhs);

    const std::vector<double> reference_x = {-2.47816,0.0797183,-0.601563,-2.47816,-3.65944,-2.77503,1.261,0.0797183,-0.0218442,-0.304688,-0.5,-0.601562};

    KRATOS_CHECK_VECTOR_NEAR(x, reference_x, 1.0e-4);
}
}
