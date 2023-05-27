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
namespace Kratos
{
namespace Testing
{

    KRATOS_TEST_CASE_IN_SUITE(HelmholtzSolidElement, KratosOptimizationFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_SCALAR);

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

        Matrix lhs;
        Vector rhs;
        VariableUtils().SetVariable(HELMHOLTZ_SCALAR_SOURCE, 1.0, r_model_part.Nodes(), ACTIVE, false);
        VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes(), ACTIVE, false);

        p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

        Vector x;
        MathUtils<double>::Solve(lhs,x,rhs);

        const std::vector<double> reference_x = {-1.0,-1.0,-1.0,-1.0};

        KRATOS_CHECK_VECTOR_NEAR(x, reference_x, 1.0e-6);
    }
    KRATOS_TEST_CASE_IN_SUITE(HelmholtzVectorSolidElement, KratosOptimizationFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart", 1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(HELMHOLTZ_VECTOR);

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
        VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes(), ACTIVE, false);
        VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes(), ACTIVE, false);

        p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

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
        VariableUtils().SetVariable(HELMHOLTZ_VECTOR_SOURCE, one_array, r_model_part.Nodes(), ACTIVE, false);
        VariableUtils().SetVariable(NUMBER_OF_NEIGHBOUR_ELEMENTS, 1.0, r_model_part.Nodes(), ACTIVE, false);

        p_element->CalculateLocalSystem(lhs,rhs,r_model_part.GetProcessInfo());

        Vector x;
        MathUtils<double>::Solve(lhs,x,rhs);

        const std::vector<double> reference_x = {-2.47816,0.0797183,-0.601563,-2.47816,-3.65944,-2.77503,1.261,0.0797183,-0.0218442,-0.304688,-0.5,-0.601562};

        KRATOS_CHECK_VECTOR_NEAR(x, reference_x, 1.0e-4);
    }
}
}
