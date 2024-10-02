// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "boost/numeric/ublas/assignment.hpp"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"
#include "custom_constitutive/plane_strain_type.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace
{

using namespace Kratos;

Vector CalculateStress(GeoLinearElasticPlaneStrain2DLaw& rConstitutiveLaw)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    auto strain = Vector{ScalarVector{4, 1.0}};
    parameters.SetStrainVector(strain);

    Vector stress;
    parameters.SetStressVector(stress);

    Matrix constitutive_matrix;
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    Properties properties;
    properties.SetValue(YOUNG_MODULUS, 1.0e7);
    properties.SetValue(POISSON_RATIO, 0.3);
    parameters.SetMaterialProperties(properties);

    rConstitutiveLaw.CalculateMaterialResponsePK2(parameters);

    return stress;
}

GeoLinearElasticPlaneStrain2DLaw CreateLinearElasticPlaneStrainLaw()
{
    return GeoLinearElasticPlaneStrain2DLaw{std::make_unique<PlaneStrainDimension>()};
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawClone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto       law   = CreateLinearElasticPlaneStrainLaw();
    const auto clone = law.Clone();
    KRATOS_EXPECT_NE(&law, clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const GeoLinearElasticPlaneStrain2DLaw*>(clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawRequiresInitializeResponse, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();
    KRATOS_EXPECT_TRUE(law.RequiresInitializeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawRequiresFinalizeResponse, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();
    KRATOS_EXPECT_TRUE(law.RequiresFinalizeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawIsIncremental, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();
    KRATOS_EXPECT_TRUE(law.IsIncremental())
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedLawFeatures, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();

    ConstitutiveLaw::Features law_features;
    law.GetLawFeatures(law_features);

    KRATOS_EXPECT_TRUE(law_features.mOptions.IsDefined(ConstitutiveLaw::PLANE_STRAIN_LAW))
    KRATOS_EXPECT_TRUE(law_features.mOptions.IsDefined(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(law_features.mOptions.IsDefined(ConstitutiveLaw::ISOTROPIC))

    const auto strain_measures = law_features.mStrainMeasures;
    KRATOS_EXPECT_TRUE(std::find(strain_measures.begin(), strain_measures.end(),
                                 ConstitutiveLaw::StrainMeasure_Infinitesimal) != strain_measures.end());
    KRATOS_EXPECT_TRUE(std::find(strain_measures.begin(), strain_measures.end(),
                                 ConstitutiveLaw::StrainMeasure_Deformation_Gradient) !=
                       strain_measures.end());

    KRATOS_EXPECT_EQ(law_features.mStrainSize, 4);
    KRATOS_EXPECT_EQ(law_features.mSpaceDimension, 2);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();
    KRATOS_EXPECT_EQ(law.GetStrainSize(), 4);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedWorkingSpaceDimension,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();
    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 2);
}

Vector CalculateStrain(GeoLinearElasticPlaneStrain2DLaw& rConstitutiveLaw)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    Vector strain;
    parameters.SetStrainVector(strain);

    rConstitutiveLaw.CalculateMaterialResponsePK2(parameters);

    return strain;
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedStress, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();

    const auto stress = CalculateStress(law);

    Vector expected_stress{4};
    expected_stress <<= 2.5e+07, 2.5e+07, 2.5e+07, 3.84615e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedStress_WhenOnlyDiagonalEntriesAreConsidered,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();
    law.SetConsiderDiagonalEntriesOnlyAndNoShear(true);

    const auto stress = CalculateStress(law);

    Vector expected_stress{4};
    expected_stress <<= 1.34615e+07, 1.34615e+07, 1.34615e+07, 0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedStress_WithInitialStressAndStrain,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();

    ConstitutiveLaw::Parameters parameters;
    auto                        initial_strain = Vector{ScalarVector{4, 0.5}};
    parameters.SetStrainVector(initial_strain);
    auto initial_stress = Vector{ScalarVector{4, 1e6}};
    parameters.SetStressVector(initial_stress);
    law.InitializeMaterialResponseCauchy(parameters);

    const auto stress = CalculateStress(law);

    Vector expected_stress{4};
    expected_stress <<= 1.35e+07, 1.35e+07, 1.35e+07, 2.92308e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedStress_AfterFinalizeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();

    ConstitutiveLaw::Parameters initial_parameters;
    auto                        initial_strain = Vector{ScalarVector{4, 0.5}};
    initial_parameters.SetStrainVector(initial_strain);
    auto initial_stress = Vector{ScalarVector{4, 1e6}};
    initial_parameters.SetStressVector(initial_stress);
    law.InitializeMaterialResponseCauchy(initial_parameters);

    auto stress = CalculateStress(law);

    ConstitutiveLaw::Parameters final_parameters;
    auto                        final_strain = Vector{ScalarVector{4, 1.3}};
    final_parameters.SetStrainVector(final_strain);
    law.FinalizeMaterialResponseCauchy(final_parameters);
    stress = CalculateStress(law);

    Vector expected_stress{4};
    expected_stress <<= 6e+06, 6e+06, 6e+06, 1.76923e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedStrain, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto                        law = CreateLinearElasticPlaneStrainLaw();
    ConstitutiveLaw::Parameters parameters;
    auto strain = Vector{3}; // The law expects the vector to have the correct length (which somehow is 3)
    parameters.SetStrainVector(strain);

    // clang-format off
    auto deformation_gradient = Matrix{3,3};
    deformation_gradient <<= 1.0, 2.0, 3.0,
                             4.0, 5.0, 6.0,
                             7.0, 8.0, 9.0;
    // clang-format on
    parameters.SetDeformationGradientF(deformation_gradient);

    law.CalculateMaterialResponsePK2(parameters);

    Vector expected_strain{3};
    expected_strain <<= 8.0, 14.0, 22.0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_strain, strain, 1e-3);
}

} // namespace Kratos::Testing