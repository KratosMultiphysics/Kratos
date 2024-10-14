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

#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"
#include "custom_constitutive/plane_strain.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

Vector CalculateStress(GeoLinearElasticPlaneStrain2DLaw& rConstitutiveLaw)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

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
    return GeoLinearElasticPlaneStrain2DLaw{std::make_unique<PlaneStrain>()};
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsCloneOfCorrectType, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law         = CreateLinearElasticPlaneStrainLaw();
    const auto p_law_clone = law.Clone();
    KRATOS_EXPECT_NE(&law, p_law_clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const GeoLinearElasticPlaneStrain2DLaw*>(p_law_clone.get()), nullptr);
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

    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::PLANE_STRAIN_LAW))
    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::ISOTROPIC))

    const auto& strain_measures = law_features.mStrainMeasures;
    KRATOS_EXPECT_NE(std::find(strain_measures.begin(), strain_measures.end(), ConstitutiveLaw::StrainMeasure_Infinitesimal),
                     strain_measures.end());
    KRATOS_EXPECT_NE(std::find(strain_measures.begin(), strain_measures.end(),
                               ConstitutiveLaw::StrainMeasure_Deformation_Gradient),
                     strain_measures.end());

    KRATOS_EXPECT_EQ(law_features.mStrainSize, 4);
    KRATOS_EXPECT_EQ(law_features.mSpaceDimension, 2);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law = CreateLinearElasticPlaneStrainLaw();
    KRATOS_EXPECT_EQ(law.GetStrainSize(), 4);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedWorkingSpaceDimension,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();
    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 2);
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

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawThrows_WhenElementProvidedStrainIsSetToFalse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateLinearElasticPlaneStrainLaw();

    ConstitutiveLaw::Parameters parameters;
    parameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        law.CalculateMaterialResponsePK2(parameters),
        "The GeoLinearElasticLaw needs an element provided strain");
}

} // namespace Kratos::Testing