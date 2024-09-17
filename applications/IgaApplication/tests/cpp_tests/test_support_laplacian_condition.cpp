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
#include "custom_conditions/support_laplacian_condition.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos::Testing
{
namespace
{
    typedef std::size_t SizeType;

    typename Condition::Pointer GetSupportLaplacianCondition(
        ModelPart& rModelPart, 
        SizeType PolynomialDegree, 
        IntegrationPoint<3> IntegrationPoint)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(PENALTY_FACTOR, -1); // Penalty-free formulation

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometryOnCurve(
            rModelPart, PolynomialDegree, IntegrationPoint);
        p_quadrature_point->SetValue(TEMPERATURE, 0.24);

        return Kratos::make_intrusive<SupportLaplacianCondition>(1, p_quadrature_point, p_elem_prop);
    }


    // Tests the stiffness matrix of the SupportLaplacianCondition with a polynomial degree of p=3.
    KRATOS_TEST_CASE_IN_SUITE(SupportLaplacianConditionP2, KratosIgaFastSuite)
    {
        Model model;
        auto &r_model_part = model.CreateModelPart("ModelPart");

        // Set buffer size
        r_model_part.SetBufferSize(2);

        // Set convection diffusion settings
        auto p_conv_dff_set = Kratos::make_shared<ConvectionDiffusionSettings>();
        p_conv_dff_set->SetDiffusionVariable(CONDUCTIVITY);
        p_conv_dff_set->SetUnknownVariable(TEMPERATURE);
        p_conv_dff_set->SetVolumeSourceVariable(HEAT_FLUX);
        r_model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_dff_set);
        // Variables addition
        r_model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);
        r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
        r_model_part.AddNodalSolutionStepVariable(HEAT_FLUX);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
        auto p_support_condition = GetSupportLaplacianCondition(r_model_part, 3, integration_point);

        for (auto& r_node : r_model_part.Nodes()) {
            r_node.AddDof(TEMPERATURE);
        }
        p_support_condition->Initialize(r_process_info);
        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

        //Check RHS and LHS results
        const double tolerance = 1.0e-6;

        const std::array<double, 8> expected_LHS{0,0,0,-2.1684e-19,0.00763467,0.011452,0.00572601,0.000954334};
        const std::array<double, 8> expected_RHS{0.00618409,0.00927613,0.00463806,0.000773011,-0.00618409,-0.00927613,-0.00463806,-0.000773011};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
        KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
        KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }
}
}
