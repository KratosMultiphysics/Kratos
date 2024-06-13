// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "testing/testing.h"

// Application includes
#include "constitutive_laws_application_variables.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_simo_taylor_neo_hookean_3d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/finite_strains/hyperelasticity/hyper_elastic_simo_taylor_neo_hookean_plane_strain_2d.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
namespace Testing
{

namespace
{
    void AuxiliaryHyperElasticConstitutiveLawTest(
        const Matrix& rDeformationGradient,
        const double YoungModulus,
        const double PoissonRatio,
        ConstitutiveLaw::Pointer pConstitutiveLaw,
        Vector& rStress,
        Matrix& rConstitutiveMatrix)
    {
        // Calculate the input strain from the deformation field
        const std::size_t strain_size = pConstitutiveLaw->GetStrainSize();
        Vector E_vector(strain_size);
        const double det_F = MathUtils<double>::Det(rDeformationGradient);
        Matrix C_matrix = prod(trans(rDeformationGradient), rDeformationGradient);
        if (strain_size == 3) {
            ConstitutiveLawUtilities<3>::CalculateGreenLagrangianStrain(C_matrix, E_vector);
        } else {
            ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_matrix, E_vector);
        }

        // Set the material properties container
        Properties material_properties;
        material_properties.SetValue(YOUNG_MODULUS, YoungModulus);
        material_properties.SetValue(POISSON_RATIO, PoissonRatio);

        // Set the constitutive law options
        Flags cl_options;
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        // Set the constitutive law input parameters
        ConstitutiveLaw::Parameters cl_parameters;
        cl_parameters.SetOptions(cl_options);
        cl_parameters.SetStrainVector(E_vector);
        cl_parameters.SetDeterminantF(det_F);
        cl_parameters.SetDeformationGradientF(rDeformationGradient);
        cl_parameters.SetMaterialProperties(material_properties);

        // Set the constitutive law output containers
        cl_parameters.SetStressVector(rStress);
        cl_parameters.SetConstitutiveMatrix(rConstitutiveMatrix);

        // Call the CL to calculate the material response
        pConstitutiveLaw->CalculateMaterialResponsePK2(cl_parameters);
    }
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawHyperElasticIsotropicNeoHookeanPlaneStrain2D, KratosConstitutiveLawsFastSuite)
{
    // Set the deformation field and test settings
    Matrix def_grad = ZeroMatrix(2,2);
    def_grad(0,0) = 2.0;
    def_grad(1,1) = 3.0/4.0;
    def_grad(1,0) = 1.0;
    const double poisson_ratio = 0.25;
    const double young_modulus = 250.0;
    auto p_neo_hookean_cl = Kratos::make_shared<HyperElasticIsotropicNeoHookeanPlaneStrain2D>();

    // Perform the material response calculation
    Matrix C_mat(3,3);
    Vector pk2_vect(3);
    AuxiliaryHyperElasticConstitutiveLawTest(def_grad, young_modulus, poisson_ratio, p_neo_hookean_cl, pk2_vect, C_mat);

    // Check the results
    std::vector<double> pk2_stress_result({85.1366277027,-32.1188648649,19.8178297297});
    std::vector<double> const_mat_row_0_result({13.6816861486,68.767442042,-18.2422481982});
    KRATOS_EXPECT_VECTOR_NEAR(pk2_vect, pk2_stress_result, 1.0e-8);
    KRATOS_EXPECT_VECTOR_NEAR(row(C_mat, 0), const_mat_row_0_result, 1.0e-8);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawHyperElasticIsotropicNeoHookean3D, KratosConstitutiveLawsFastSuite)
{
    // Set the deformation field and test settings
    Matrix def_grad = ZeroMatrix(3,3);
    def_grad(0,0) = 2.0;
    def_grad(1,1) = 3.0/4.0;
    def_grad(2,2) = 0.5;
    def_grad(0,2) = 0.5;
    def_grad(1,0) = 1.0;
    def_grad(2,0) = 1.0;
    const double poisson_ratio = 0.25;
    const double young_modulus = 250.0;
    auto p_neo_hookean_cl = Kratos::make_shared<HyperElasticIsotropicNeoHookean3D>();

    // Perform the material response calculation
    Matrix C_mat(6,6);
    Vector pk2_vect(6);
    AuxiliaryHyperElasticConstitutiveLawTest(def_grad, young_modulus, poisson_ratio, p_neo_hookean_cl, pk2_vect, C_mat);

    // Check the results
    std::vector<double> pk2_stress_result({-296.16585060235,-956.44226827292,-3861.6585060235,528.22113413646,-1584.6634024094,1188.497551807});
    std::vector<double> const_mat_row_0_result({1984.6634024094,3883.8460487278,18261.970621684,-2646.2178698792,7938.6536096375,-5953.9902072281});
    KRATOS_EXPECT_VECTOR_NEAR(pk2_vect, pk2_stress_result, 1.0e-8);
    KRATOS_EXPECT_VECTOR_NEAR(row(C_mat, 0), const_mat_row_0_result, 1.0e-8);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawHyperElasticSimoTaylorNeoHookeanPlaneStrain2D, KratosConstitutiveLawsFastSuite)
{
    // Set the deformation field and test settings
    Matrix def_grad = ZeroMatrix(2,2);
    def_grad(0,0) = 2.0;
    def_grad(1,1) = 3.0/4.0;
    def_grad(1,0) = 1.0;
    const double young_modulus = 250.0;
    const double poisson_ratio = 0.49995;
    auto p_simo_taylor_cl = Kratos::make_shared<HyperElasticSimoTaylorNeoHookeanPlaneStrain2D>();

    // Perform the material response calculation
    Matrix C_mat(3,3);
    Vector pk2_vect(3);
    AuxiliaryHyperElasticConstitutiveLawTest(def_grad, young_modulus, poisson_ratio, p_simo_taylor_cl, pk2_vect, C_mat);

    // Check the results
    std::vector<double> pk2_stress_result({130225.260981,1157119.58917,-173559.604765});
    std::vector<double> const_mat_row_0_result({52084.5269495,925908.757144,-69464.5550685});
    KRATOS_EXPECT_VECTOR_NEAR(pk2_vect, pk2_stress_result, 1.0e-5);
    KRATOS_EXPECT_VECTOR_NEAR(row(C_mat, 0), const_mat_row_0_result, 1.0e-5);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawHyperElasticSimoTaylorNeoHookean3D, KratosConstitutiveLawsFastSuite)
{
    // Set the deformation field and test settings
    Matrix def_grad = ZeroMatrix(3,3);
    def_grad(0,0) = 2.0;
    def_grad(1,1) = 3.0/4.0;
    def_grad(2,2) = 0.5;
    def_grad(0,2) = 0.5;
    def_grad(1,0) = 1.0;
    def_grad(2,0) = 1.0;
    const double young_modulus = 250.0;
    const double poisson_ratio = 0.49995;
    auto p_simo_taylor_cl = Kratos::make_shared<HyperElasticSimoTaylorNeoHookean3D>();

    // Perform the material response calculation
    Matrix C_mat(6,6);
    Vector pk2_vect(6);
    AuxiliaryHyperElasticConstitutiveLawTest(def_grad, young_modulus, poisson_ratio, p_simo_taylor_cl, pk2_vect, C_mat);

    // Check the results
    std::vector<double> pk2_stress_result({-716740.113382,-1911574.06118,-7168843.43147,955867.158235,-2867601.47471,2150701.10603});
    std::vector<double> const_mat_row_0_result({3336930.174379,6349857.4992591,30503623.332614,-4449525.130806,13348575.392418,-10011431.544314});
    KRATOS_EXPECT_VECTOR_NEAR(pk2_vect, pk2_stress_result, 1.0e-5);
    KRATOS_EXPECT_VECTOR_NEAR(row(C_mat, 0), const_mat_row_0_result, 1.0e-5);
}

} // namespace Testing
} // namespace Kratos
