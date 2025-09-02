// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#include "custom_constitutive/small_strain_umat_2D_plane_strain_law.hpp"
#include "custom_constitutive/three_dimensional.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <string>

using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(SmallStrainUMATLaw_HasReturnsTrueForCauchyStress,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto umat_law = SmallStrainUMAT2DPlaneStrainLaw{std::make_unique<ThreeDimensional>()};

    // Act & Assert
    KRATOS_EXPECT_TRUE(umat_law.Has(STATE_VARIABLES));
    KRATOS_EXPECT_TRUE(umat_law.Has(CAUCHY_STRESS_VECTOR));
    KRATOS_EXPECT_FALSE(umat_law.Has(TOTAL_STRESS_VECTOR));
}



} // namespace Kratos::Testing