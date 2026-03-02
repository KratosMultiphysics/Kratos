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

#include "custom_constitutive/small_strain_umat_2D_plane_strain_law.h"
#include "custom_constitutive/three_dimensional.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

class ParametrizedUMATLawTests : public ::testing::TestWithParam<std::tuple<Variable<Vector>, bool>>
{
};

TEST_P(ParametrizedUMATLawTests, SmallStrainUMATLaw_HasReturnsCorrectBoolValueForVectorVariables)
{
    // Arrange
    const auto& [vector_variable, expected_state] = GetParam();
    auto umat_law = SmallStrainUMAT2DPlaneStrainLaw{std::make_unique<ThreeDimensional>()};

    // Act & Assert
    KRATOS_EXPECT_EQ(umat_law.Has(vector_variable), expected_state);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedUMATLawTests,
                         ::testing::Values(std::make_tuple(STATE_VARIABLES, true),
                                           std::make_tuple(CAUCHY_STRESS_VECTOR, true),
                                           std::make_tuple(TOTAL_STRESS_VECTOR, false)));
} // namespace Kratos::Testing