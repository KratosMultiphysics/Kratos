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
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometry(
            rModelPart, PolynomialDegree, IntegrationPoint);

        return Kratos::make_intrusive<Shell3pElement>(1, p_quadrature_point, p_elem_prop);
    }


    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=2.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP2, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0, 0.0, 0.0, 1.0);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 2, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-5;
        const std::vector<double> expected_right_hand_side_vector({7.83999e-08,2.07838e-11,2.85724e-10,1.19805e-07,2.53715e-09,6.58355e-10,6.34008e-08,5.77027e-11,3.9188e-10,7.30737e-08,2.26601e-10,2.88665e-10,1.56401e-10,6.67553e-08,2.00164e-10,3.7924e-09,4.32224e-07,2.78981e-09,3.98132e-10,5.89574e-08,2.76404e-10,1.40101e-10,3.35175e-08,1.16888e-10,5.6946e-10,5.4851e-10,1.3056e-07,2.87022e-10,2.06904e-12,3.2299e-05,4.79543e-10,9.35669e-11,1.21407e-07,1.03257e-09,2.52929e-10,6.84278e-08,1.20185e-07,1.68722e-08,9.15427e-10,1.88326e-07,3.92703e-10,5.80919e-11,3.79004e-07,7.25958e-10,8.25508e-10,2.00979e-07,1.47245e-09,1.6338e-08,1.42025e-09,4.38993e-07,6.80859e-09,2.92253e-11,1.19192e-07,4.58938e-11,1.20997e-10,7.48042e-08,1.48577e-10,3.63387e-09,7.41823e-08,1.16976e-11,1.08821e-09,8.9436e-09,3.25153e-05,3.11555e-10,4.59753e-11,3.59894e-07,2.29186e-09,6.19311e-10,2.72257e-07,1.00945e-08,1.85752e-09,2.25757e-07,6.36167e-08,2.45672e-10,9.80188e-10,3.78806e-07,1.66561e-10,3.77833e-10,1.51191e-07,1.15385e-10,4.436e-10,8.086e-06,4.61353e-11,1.95153e-10,2.45672e-10,5.88961e-08,1.05936e-09,4.80807e-11,7.45737e-08,1.67813e-10,2.08987e-10,1.12515e-07,5.10433e-11,1.02384e-09,1.9622e-07,7.09616e-10,1.02271e-10,1.0528e-09,1.21355e-07,7.72253e-10,3.54967e-10,2.72712e-07,4.1818e-10,1.81686e-10,3.05914e-07,8.06128e-10,1.81512e-09,3.57271e-07,7.28385e-08,7.78619e-10,3.12537e-10,2.00335e-07,3.25997e-09,4.32728e-09,8.16507e-06,1.47387e-09,7.04346e-10,1.60776e-07,1.28813e-09,1.48714e-09,3.26769e-10,3.37628e-08,6.24457e-10,7.40437e-10,7.46591e-08,9.46029e-11,3.09336e-10,1.9826e-07,9.3412e-10,1.80188e-10,7.87286e-08,2.51837e-10,1.04632e-09,4.97367e-10,6.78623e-08,7.14672e-09,3.06683e-10,2.24474e-07,2.04464e-09,3.84428e-09,3.56122e-07,1.10454e-09,3.60626e-10,1.84727e-07});
        KRATOS_CHECK_MATRIX_NEAR(left_hand_side_matrix, left_hand_side_matrix, tolerance);
        KRATOS_CHECK_VECTOR_NEAR(right_hand_side_vector, expected_right_hand_side_vector, tolerance);
    }
}
}
