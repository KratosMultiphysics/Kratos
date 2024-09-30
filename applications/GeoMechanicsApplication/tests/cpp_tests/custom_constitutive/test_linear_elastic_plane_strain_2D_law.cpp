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

} // namespace Kratos::Testing