//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "test_creation_utility.h"
#include "custom_conditions/support_pressure_condition.h"

namespace Kratos::Testing
{
namespace
{
    typedef std::size_t SizeType;

    typename Condition::Pointer GetSupportPressureCondition(
        ModelPart& rModelPart, 
        SizeType PolynomialDegree, 
        IntegrationPoint<3> IntegrationPoint)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometryOnCurve(
            rModelPart, PolynomialDegree, IntegrationPoint);
        p_quadrature_point->SetValue(PRESSURE, 2500.0);
        Vector normal_stress = ZeroVector(2);
        normal_stress[0] = 1000.0; // Normal stress in x-direction
        normal_stress[1] = -1000.0; // Normal stress in y-direction
        p_quadrature_point->SetValue(NORMAL_STRESS, normal_stress);

        Vector mesh_size(2); 
        mesh_size[0] = 0.1;  
        mesh_size[1] = 0.1;
        p_quadrature_point->SetValue(KNOT_SPAN_SIZES, mesh_size);

        return Kratos::make_intrusive<SupportPressureCondition>(1, p_quadrature_point, p_elem_prop);
    }
}


// Tests the stiffness matrix of the SupportLaplacianCondition with a polynomial degree of p=3.
KRATOS_TEST_CASE_IN_SUITE(SupportPressureConditionP2, KratosIgaFastSuite)
{
    Model model;
    auto &r_model_part = model.CreateModelPart("ModelPart");

    // Set buffer size
    r_model_part.SetBufferSize(2);

    // Variables addition
    r_model_part.AddNodalSolutionStepVariable(VELOCITY_X);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY_Y);
    r_model_part.AddNodalSolutionStepVariable(PRESSURE);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

    const auto& r_process_info = r_model_part.GetProcessInfo();

    IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
    auto p_support_condition = GetSupportPressureCondition(r_model_part, 3, integration_point);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
        r_node.AddDof(PRESSURE);
    }
    p_support_condition->Initialize(r_process_info);
    Matrix left_hand_side_matrix;
    Vector right_hand_side_vector;
    p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

    //Check RHS and LHS results
    const double tolerance = 1.0e-7;
    const std::array<double, 24> expected_LHS{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    const std::array<double, 24> expected_RHS{24.4786743,36.7180114,0.0000000,36.7180114,55.0770171,0.0000000,18.3590057,27.5385086,0.0000000,3.0598343,4.5897514,0.0000000,1.2883513,1.9325269,0.0000000,1.9325269,2.8987904,0.0000000,0.9662635,1.4493952,0.0000000,0.1610439,0.2415659,0.0000000};

    for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
    KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
    }
    for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
    KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
    }
}
}
