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
//                   Gennady Markelov
//

#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/three_dimensional.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

namespace
{

using namespace Kratos;

Vector Calculate3DStress(GeoIncrementalLinearElasticLaw& rConstitutiveLaw)
{
    ConstitutiveLaw::Parameters parameters;
    parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);

    auto strain = Vector{ScalarVector{6, 1.0}};
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

GeoIncrementalLinearElasticLaw CreateIncrementalLinearElastic3DLaw()
{
    return GeoIncrementalLinearElasticLaw{std::make_unique<ThreeDimensional>()};
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedLawFeatures, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();

    ConstitutiveLaw::Features law_features;
    law.GetLawFeatures(law_features);

    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::THREE_DIMENSIONAL_LAW))
    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::INFINITESIMAL_STRAINS))
    KRATOS_EXPECT_TRUE(law_features.mOptions.Is(ConstitutiveLaw::ISOTROPIC))

    const auto& strain_measures = law_features.mStrainMeasures;
    KRATOS_EXPECT_NE(std::find(strain_measures.begin(), strain_measures.end(), ConstitutiveLaw::StrainMeasure_Infinitesimal),
                     strain_measures.end());
    KRATOS_EXPECT_NE(std::find(strain_measures.begin(), strain_measures.end(),
                               ConstitutiveLaw::StrainMeasure_Deformation_Gradient),
                     strain_measures.end());

    KRATOS_EXPECT_EQ(law_features.mStrainSize, 6);
    KRATOS_EXPECT_EQ(law_features.mSpaceDimension, 3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedStrainSize, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law = CreateIncrementalLinearElastic3DLaw();
    KRATOS_EXPECT_EQ(law.GetStrainSize(), 6);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedWorkingSpaceDimension,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();
    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedStress, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();

    const auto stress = Calculate3DStress(law);

    Vector expected_stress{6};
    expected_stress <<= 2.5e+07, 2.5e+07, 2.5e+07, 3.84615e+06, 3.84615e+06, 3.84615e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedStress_WhenOnlyDiagonalEntriesAreConsidered,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();
    law.SetConsiderDiagonalEntriesOnlyAndNoShear(true);

    const auto stress = Calculate3DStress(law);

    Vector expected_stress{6};
    expected_stress <<= 1.34615e+07, 1.34615e+07, 1.34615e+07, 0, 0, 0;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedStress_WithInitialStressAndStrain,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();

    ConstitutiveLaw::Parameters parameters;
    auto                        initial_strain = Vector{ScalarVector{6, 0.5}};
    parameters.SetStrainVector(initial_strain);
    auto initial_stress = Vector{ScalarVector{6, 1e6}};
    parameters.SetStressVector(initial_stress);
    law.InitializeMaterialResponseCauchy(parameters);

    const auto stress = Calculate3DStress(law);

    Vector expected_stress{6};
    expected_stress <<= 1.35e+07, 1.35e+07, 1.35e+07, 2.92308e+06, 2.92308e+06, 2.92308e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedStress_AfterFinalizeMaterialResponse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();

    ConstitutiveLaw::Parameters initial_parameters;
    auto                        initial_strain = Vector{ScalarVector{6, 0.5}};
    initial_parameters.SetStrainVector(initial_strain);
    auto initial_stress = Vector{ScalarVector{6, 1e6}};
    initial_parameters.SetStressVector(initial_stress);
    law.InitializeMaterialResponseCauchy(initial_parameters);

    auto stress = Calculate3DStress(law);

    ConstitutiveLaw::Parameters final_parameters;
    auto                        final_strain = Vector{ScalarVector{6, 1.3}};
    final_parameters.SetStrainVector(final_strain);
    law.FinalizeMaterialResponseCauchy(final_parameters);
    stress = Calculate3DStress(law);

    Vector expected_stress{6};
    expected_stress <<= 6e+06, 6e+06, 6e+06, 1.76923e+06, 1.76923e+06, 1.76923e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawReturnsExpectedStress_AfterFinalizeMaterialResponseAndReset,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();

    ConstitutiveLaw::Parameters initial_parameters;
    auto                        initial_strain = Vector{ ScalarVector{6, 0.5} };
    initial_parameters.SetStrainVector(initial_strain);
    auto initial_stress = Vector{ ScalarVector{6, 1e6} };
    initial_parameters.SetStressVector(initial_stress);
    law.InitializeMaterialResponseCauchy(initial_parameters);

    auto stress = Calculate3DStress(law);

    ConstitutiveLaw::Parameters final_parameters;
    auto                        final_strain = Vector{ ScalarVector{6, 1.3} };
    final_parameters.SetStrainVector(final_strain);
    law.FinalizeMaterialResponseCauchy(final_parameters);
    
    const Properties properties;
    const Geometry<Node> geometry;
    const Vector shape_functions_values;

    law.ResetMaterial(properties, geometry, shape_functions_values);

    stress = Calculate3DStress(law);

    Vector expected_stress{ 6 };
    expected_stress <<= 2.5e+07, 2.5e+07, 2.5e+07, 3.84615e+06, 3.84615e+06, 3.84615e+06;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress, stress, 1e-3);
}

#ifdef KRATOS_DEBUG
KRATOS_TEST_CASE_IN_SUITE(GeoIncrementalLinearElastic3DLawThrows_WhenElementProvidedStrainIsSetToFalse,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto law = CreateIncrementalLinearElastic3DLaw();

    ConstitutiveLaw::Parameters parameters;
    parameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(law.CalculateMaterialResponsePK2(parameters),
                                      "The GeoLinearElasticLaw needs an element provided strain");
}
#endif

} // namespace Kratos::Testing