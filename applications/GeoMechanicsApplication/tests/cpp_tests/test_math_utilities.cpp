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

#include "custom_utilities/math_utilities.h"
#include "includes/checks.h"
#include "geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CalculateDeterminants_ReturnsCorrectResults, KratosGeoMechanicsFastSuite)
{
    const std::vector<Matrix> matrices = {ScalarMatrix(1, 1, 3.0), ZeroMatrix(3, 3), IdentityMatrix(3) * 2.0};

    const std::vector<double> results = GeoMechanicsMathUtilities::CalculateDeterminants(matrices);

    const std::vector<double> expected_results = {3.0, 0.0, 8.0};
    KRATOS_CHECK_VECTOR_NEAR(results, expected_results, 1.0e-12)
}

KRATOS_TEST_CASE_IN_SUITE(CalculateDeterminants_ReturnsEmptyVectorForEmptyInput, KratosGeoMechanicsFastSuite)
{
    const std::vector<Matrix> matrices = {};

    const std::vector<double> results = GeoMechanicsMathUtilities::CalculateDeterminants(matrices);

    KRATOS_EXPECT_TRUE(results.empty())
}

} // namespace Kratos::Testing