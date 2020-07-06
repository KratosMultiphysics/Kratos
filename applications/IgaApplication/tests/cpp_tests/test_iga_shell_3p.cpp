//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

#include "test_creation_utility.h"

#include "custom_elements/shell_3p_element.h"

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

    typename Shell3pElement::Pointer GetShell3pElement(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 200000000);
        p_elem_prop->SetValue(POISSON_RATIO, 0.0);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometry(
            rModelPart, PolynomialDegree, IntegrationPoint);

        return Kratos::make_intrusive<Shell3pElement>(1, p_quadrature_point, p_elem_prop);
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=3.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP3, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0694318, 0.211324865405187, 0.0, 0.0869637);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 3, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-1;

        const std::vector<double> expected_LHS_row_0({637725.123769309,143581.472497936,0,64285.8786450023,-122155.580440902,0,-1057.60612106728,-20626.5743552912,0,-171.900069634590,-799.317701742961,0,-545148.472653437,38472.5396038933,0,-143047.813284562,-32731.4891300946,0,-12241.7363677962,-5526.87394112079,0,-343.473917813981,-214.176532677927,0});
        const std::vector<double> expected_right_hand_side_vector({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
        
        KRATOS_CHECK_VECTOR_NEAR(row(left_hand_side_matrix,0), expected_LHS_row_0, tolerance);
        KRATOS_CHECK_VECTOR_NEAR(right_hand_side_vector, expected_right_hand_side_vector, tolerance);
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=4.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP4, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0469100770306680, 0.211324865405187, 0, 0.0592317212640473);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 4, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-1;

        const std::vector<double> expected_LHS_row_0({491667.161414150,133490.278045005,0,4078.07575802542,-113779.527542805,0,-6544.23883310685,-18740.6082290450,0,-439.346498780415,-954.225814603370,0,-8.16984378284357,-15.9164585515075,0,-379618.738844226,35768.6121995651,0,-99581.7841076655,-30487.1325202893,0,-9186.16206812944,-5021.53084064066,0,-361.606940433798,-255.684036419903,0,-5.19003605046336,-4.26480221523988,0});
        const std::vector<double> expected_right_hand_side_vector({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
        
        KRATOS_CHECK_VECTOR_NEAR(row(left_hand_side_matrix,0), expected_LHS_row_0, tolerance);
        KRATOS_CHECK_VECTOR_NEAR(right_hand_side_vector, expected_right_hand_side_vector, tolerance);
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=5.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP5, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0337652428984240, 0.211324865405187, 0, 0.0428311230947926);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 5, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-1;

        //Matrix expected_left_hand_side_matrix = ZeroMatrix(36,36);
        const std::vector<double> expected_LHS_row_0({405000.290031407,123985.676772767,0,-33973.9881469315,-106654.871970357,0,-9694.63076821954,-16422.3627864133,0,-594.585433000862,-887.278239041373,0,-14.8585189361368,-20.9788839821913,0,-0.135084710956110,-0.184892973551875,0,-276681.912058038,33221.8619642891,0,-76407.9570133629,-28578.0868133018,0,-7301.61436342488,-4400.35884643036,0,-323.698899579919,-237.745487612845,0,-6.85347582290271,-5.62127502113436,0,-0.0562693799549848,-0.0495419229494137,0});
        const std::vector<double> expected_right_hand_side_vector({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0});
        
        KRATOS_CHECK_VECTOR_NEAR(row(left_hand_side_matrix,0), expected_LHS_row_0, tolerance);
        KRATOS_CHECK_VECTOR_NEAR(right_hand_side_vector, expected_right_hand_side_vector, tolerance);
    }
}
}
