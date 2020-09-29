// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

// Application includes
#include "custom_constitutive/user_provided_linear_elastic_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

namespace Testing
{
    /**
    * Checks the user provided linear elastic law in 2D
    */
    KRATOS_TEST_CASE_IN_SUITE(UserProvidedLinearElasticLaw2D, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        // Set the element properties
        auto p_prop = r_model_part.CreateNewProperties(0);
        p_prop->SetValue(ELASTICITY_TENSOR, IdentityMatrix(3,3));

        // Create the constitutive law
        UserProvidedLinearElasticLaw<2> user_provided_claw = UserProvidedLinearElasticLaw<2>();

        // Create the test geometry
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        Triangle2D3<Node<3>> geometry = Triangle2D3<Node<3>>(p_node_1, p_node_2, p_node_3);

        // Compute RHS and LHS
        Vector strain_vect(3);
        strain_vect[0] = 1.0;
        strain_vect[1] = 1.0;
        strain_vect[2] = 1.0;
        Vector stress_vect(3);
        Matrix constitutive_tensor;

        Flags cl_options;
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        ConstitutiveLaw::Parameters cl_parameters;
        cl_parameters.SetElementGeometry(geometry);
        cl_parameters.SetMaterialProperties(*p_prop);
        cl_parameters.SetStrainVector(strain_vect);
        cl_parameters.SetStressVector(stress_vect);
        cl_parameters.SetConstitutiveMatrix(constitutive_tensor);
        cl_parameters.SetOptions(cl_options);

        user_provided_claw.CalculateMaterialResponseCauchy(cl_parameters);

        // Check results
        const double tolerance = 1.0e-12;
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(cl_parameters.GetStressVector(), strain_vect, tolerance)
    }

    /**
    * Checks the user provided linear elastic law in 3D
    */
    KRATOS_TEST_CASE_IN_SUITE(UserProvidedLinearElasticLaw3D, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);

        // Set the element properties
        auto p_prop = r_model_part.CreateNewProperties(0);
        p_prop->SetValue(ELASTICITY_TENSOR, IdentityMatrix(6,6));

        // Create the constitutive law
        UserProvidedLinearElasticLaw<3> user_provided_claw = UserProvidedLinearElasticLaw<3>();

        // Create the test geometry
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
        auto p_node_3 = r_model_part.CreateNewNode(3, 0.0 , 1.0 , 0.0);
        auto p_node_4 = r_model_part.CreateNewNode(4, 0.0 , 0.0 , 1.0);
        Tetrahedra3D4<Node<3>> geometry = Tetrahedra3D4<Node<3>>(p_node_1, p_node_2, p_node_3, p_node_4);

        // Compute RHS and LHS
        Vector strain_vect(6);
        strain_vect[0] = 1.0;
        strain_vect[1] = 1.0;
        strain_vect[2] = 1.0;
        strain_vect[3] = 1.0;
        strain_vect[4] = 1.0;
        strain_vect[5] = 1.0;
        Vector stress_vect(6);
        Matrix constitutive_tensor;

        Flags cl_options;
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        ConstitutiveLaw::Parameters cl_parameters;
        cl_parameters.SetElementGeometry(geometry);
        cl_parameters.SetMaterialProperties(*p_prop);
        cl_parameters.SetStrainVector(strain_vect);
        cl_parameters.SetStressVector(stress_vect);
        cl_parameters.SetConstitutiveMatrix(constitutive_tensor);
        cl_parameters.SetOptions(cl_options);

        user_provided_claw.CalculateMaterialResponseCauchy(cl_parameters);

        // Check results
        const double tolerance = 1.0e-12;
        KRATOS_CHECK_VECTOR_RELATIVE_NEAR(cl_parameters.GetStressVector(), strain_vect, tolerance)
    }

} // namespace Testing
} // namespace Kratos.
