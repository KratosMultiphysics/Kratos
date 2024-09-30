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
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawClone, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoLinearElasticPlaneStrain2DLaw law;
    const auto                       clone = law.Clone();
    KRATOS_EXPECT_NE(&law, clone.get());
    KRATOS_EXPECT_NE(dynamic_cast<const GeoLinearElasticPlaneStrain2DLaw*>(clone.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawRequiresInitializeResponse, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoLinearElasticPlaneStrain2DLaw law;
    KRATOS_EXPECT_TRUE(law.RequiresInitializeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawRequiresFinalizeResponse, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoLinearElasticPlaneStrain2DLaw law;
    KRATOS_EXPECT_TRUE(law.RequiresFinalizeMaterialResponse())
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedLawFeatures, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoLinearElasticPlaneStrain2DLaw law;

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
    GeoLinearElasticPlaneStrain2DLaw law;
    KRATOS_EXPECT_EQ(law.GetStrainSize(), 4);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLinearElasticPlaneStrain2DLawReturnsExpectedWorkingSpaceDimension, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    GeoLinearElasticPlaneStrain2DLaw law;
    KRATOS_EXPECT_EQ(law.WorkingSpaceDimension(), 2);
}

/*

public:
[x] ConstitutiveLaw::Pointer Clone() const override;
[x] bool RequiresInitializeMaterialResponse() override;
[ ] void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
[x] bool RequiresFinalizeMaterialResponse() override;
[ ] void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
[ ] void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
[x] void GetLawFeatures(Features& rFeatures) override;
[x] SizeType WorkingSpaceDimension() override;
[x] SizeType GetStrainSize() const override;
[ ] bool IsIncremental() override;
[ ] bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;


*/

} // namespace Kratos::Testing