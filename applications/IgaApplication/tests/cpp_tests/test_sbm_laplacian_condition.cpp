//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
<<<<<<< HEAD
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
=======
//  License: BSD License
//           Kratos default license: kratos/license.txt
>>>>>>> 3D_sbm
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
<<<<<<< HEAD
#include "test_creation_utility.h"
#include "custom_conditions/sbm_laplacian_condition.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos::Testing
{
namespace
{
    typedef std::size_t SizeType;
=======

#include "test_creation_utility.h"

#include "custom_conditions/sbm_laplacian_condition.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{
namespace Testing
{
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Operations
    ///@{
>>>>>>> 3D_sbm

    typename Condition::Pointer GetSBMLaplacianCondition(
          ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint, 
          double penalty, std::string boundary_condition_type)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(PENALTY_FACTOR, penalty); // Penalty-free formulation

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometryOnCurve(
            rModelPart, PolynomialDegree, IntegrationPoint);
        
        rModelPart.CreateNewNode(1000, 0.0, 0.1, 0.0);
        rModelPart.CreateNewNode(1001, 0.15, 0.2, 0.0);
        rModelPart.CreateNewNode(1002, 0.25, 1.3, 0.0);
        auto& node1 = rModelPart.GetNode(1000);
        auto& node2 = rModelPart.GetNode(1001);
        auto& node3 = rModelPart.GetNode(1002);
        node1.SetValue(TEMPERATURE, 99.0);
        node2.SetValue(TEMPERATURE, 199.0);
        node3.SetValue(TEMPERATURE, 1199.0);
<<<<<<< HEAD
        node1.SetValue(HEAT_FLUX, 1.0); // Set the HEAT_FLUX for Neumann Conditions
=======
>>>>>>> 3D_sbm
        Properties::Pointer p_cond_prop = rModelPart.pGetProperties(0);

        Condition::Pointer p_cond1 = rModelPart.CreateNewCondition("LineCondition2D2N", 1, {{1000, 1001}}, p_cond_prop );
        Condition::Pointer p_cond2 = rModelPart.CreateNewCondition("LineCondition2D2N", 2, {{1001, 1002}}, p_cond_prop );
        p_quadrature_point->SetValue(NEIGHBOUR_CONDITIONS, GlobalPointersVector<Condition>({p_cond1, p_cond2}));
        p_quadrature_point->SetValue(BOUNDARY_CONDITION_TYPE, boundary_condition_type);

        Vector mesh_size(2); 
        mesh_size[0] = 0.1;  
        mesh_size[1] = 0.1;
        p_quadrature_point->SetValue(MARKER_MESHES, mesh_size);

        return Kratos::make_intrusive<SBMLaplacianCondition>(1, p_quadrature_point, p_elem_prop);
    }
<<<<<<< HEAD
}


// Tests the stiffness matrix of the SBMLaplacianCondition
KRATOS_TEST_CASE_IN_SUITE(SBMLaplacianConditionP2_penaltyFree, KratosIgaFastSuite)
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

    double penalty = -1;

    IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
    std::string boundary_condition_type = "dirichlet";

    auto p_support_condition = GetSBMLaplacianCondition(r_model_part, 3, integration_point, penalty, boundary_condition_type);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }
    p_support_condition->Initialize(r_process_info);
    Matrix left_hand_side_matrix;
    Vector right_hand_side_vector;
    p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

    //Check RHS and LHS results
    const double tolerance = 1.0e-4;

    const std::array<double, 8> expected_LHS{-0.403943,-0.225679,0.0772782,0.044566,-0.0407768,-0.0411528,-0.0105702,-9.40019e-05};
    const std::array<double, 8> expected_RHS{-59.4368,-34.6289,9.94865,6.20196,-4.33659,-3.63508,-0.38264,0.175377};

    for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
      KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
    }
    for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
      KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
    }
}

// Tests the stiffness matrix of the SBMLaplacianCondition
KRATOS_TEST_CASE_IN_SUITE(SBMLaplacianConditionP2_penalty, KratosIgaFastSuite)
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

    double penalty = 10.0;

    IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
    std::string boundary_condition_type = "dirichlet";

    auto p_support_condition = GetSBMLaplacianCondition(r_model_part, 3, integration_point, penalty, boundary_condition_type);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }
    p_support_condition->Initialize(r_process_info);
    Matrix left_hand_side_matrix;
    Vector right_hand_side_vector;
    p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

    //Check RHS and LHS results
    const double tolerance = 1.0e-3;

    const std::array<double, 8> expected_LHS{3.77237,2.03571,-0.793567,-0.434165,0.464875,0.506637,0.15798,0.0104404};
    const std::array<double, 8> expected_RHS{566.308,304.199,-120.532,-65.5272,71.4262,78.4413,24.8716,1.75377};

    for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
      KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
    }
    for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
      KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
    }
}


// Tests the stiffness matrix of the SBMLaplacianCondition -> NEUMANN TYPE
KRATOS_TEST_CASE_IN_SUITE(SBMLaplacianConditionP2_Neumann, KratosIgaFastSuite)
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

    double penalty = -1;

    IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
    std::string boundary_condition_type = "neumann";

    auto p_support_condition = GetSBMLaplacianCondition(r_model_part, 3, integration_point, penalty, boundary_condition_type);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }
    p_support_condition->Initialize(r_process_info);
    Matrix left_hand_side_matrix;
    Vector right_hand_side_vector;
    p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

    //Check RHS and LHS results
    const double tolerance = 1.0e-4;

    const std::array<double, 8> expected_LHS{-0.0189816,-0.00803494,0.00620128,0.00273667,0.00463952,0.00803494,0.0045553,0.000848855};
    const std::array<double, 8> expected_RHS{0.0125166,0.0187749,0.00938743,0.00156457,0.000658767,0.000988151,0.000494075,8.23459e-05};

    for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
      KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
    }
    for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
      KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
    }
=======


    // Tests the stiffness matrix of the SBMLaplacianCondition
    KRATOS_TEST_CASE_IN_SUITE(SBMLaplacianConditionP2_penaltyFree, KratosIgaFastSuite)
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

        double penalty = -1;

        IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
        std::string boundary_condition_type = "dirichlet";

        auto p_support_condition = GetSBMLaplacianCondition(r_model_part, 3, integration_point, penalty, boundary_condition_type);

        for (auto& r_node : r_model_part.Nodes()) {
            r_node.AddDof(TEMPERATURE);
        }
        p_support_condition->Initialize(r_process_info);
        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

        //Check RHS and LHS results
        const double tolerance = 1.0e-4;

        const std::array<double, 8> expected_LHS{-0.403943,-0.225679,0.0772782,0.044566,-0.0407768,-0.0411528,-0.0105702,-9.40019e-05};
        const std::array<double, 8> expected_RHS{-59.4368,-34.6289,9.94865,6.20196,-4.33659,-3.63508,-0.38264,0.175377};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the SBMLaplacianCondition
    KRATOS_TEST_CASE_IN_SUITE(SBMLaplacianConditionP2_penalty, KratosIgaFastSuite)
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

        double penalty = 10.0;

        IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
        std::string boundary_condition_type = "dirichlet";

        auto p_support_condition = GetSBMLaplacianCondition(r_model_part, 3, integration_point, penalty, boundary_condition_type);

        for (auto& r_node : r_model_part.Nodes()) {
            r_node.AddDof(TEMPERATURE);
        }
        p_support_condition->Initialize(r_process_info);
        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

        //Check RHS and LHS results
        const double tolerance = 1.0e-3;

        const std::array<double, 8> expected_LHS{3.77237,2.03571,-0.793567,-0.434165,0.464875,0.506637,0.15798,0.0104404};
        const std::array<double, 8> expected_RHS{566.308,304.199,-120.532,-65.5272,71.4262,78.4413,24.8716,1.75377};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }


    // Tests the stiffness matrix of the SBMLaplacianCondition -> NEUMANN TYPE
    KRATOS_TEST_CASE_IN_SUITE(SBMLaplacianConditionP2_Neumann, KratosIgaFastSuite)
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

        double penalty = -1;

        IntegrationPoint<3> integration_point(0.333333333333333, 0.05, 0.0, 0.086963711284364);
        std::string boundary_condition_type = "neumann";

        auto p_support_condition = GetSBMLaplacianCondition(r_model_part, 3, integration_point, penalty, boundary_condition_type);

        for (auto& r_node : r_model_part.Nodes()) {
            r_node.AddDof(TEMPERATURE);
        }
        p_support_condition->Initialize(r_process_info);
        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_support_condition->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_process_info);

        //Check RHS and LHS results
        const double tolerance = 1.0e-4;

        const std::array<double, 8> expected_LHS{-0.0189816,-0.00803494,0.00620128,0.00273667,0.00463952,0.00803494,0.0045553,0.000848855};
        const std::array<double, 8> expected_RHS{0.00107745,0.00161618,0.000808089,0.000134682,5.6708e-05,8.5062e-05,4.2531e-05,7.0885e-06};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    
>>>>>>> 3D_sbm
}
}
