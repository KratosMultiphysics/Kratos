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

/*

public:
[x] ConstitutiveLaw::Pointer Clone() const override;
[x] bool RequiresInitializeMaterialResponse() override;
[ ] void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
[x] bool RequiresFinalizeMaterialResponse() override;
[ ] void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
[ ] void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;
[ ] void GetLawFeatures(Features& rFeatures) override;
[ ] SizeType WorkingSpaceDimension() override;
[ ] SizeType GetStrainSize() const override;
[ ] bool IsIncremental() override;
[ ] bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;

protected:
[ ] void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;
[ ] void CalculatePK2Stress(const Vector&                rStrainVector,
[ ] Vector&                      rStressVector,
[ ] ConstitutiveLaw::Parameters& rValues) override;

*/

} // namespace Kratos::Testing