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

#include "testing/testing.h"
#include "custom_utilities/scheme_factory.hpp"
#include "spaces/ublas_space.h"

using namespace Kratos;

namespace Kratos::Testing
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using SchemeFactoryType = SchemeFactory<SparseSpaceType, LocalSpaceType>;

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_Throws_WhenSchemeTypeIsUndefined, KratosGeoMechanicsFastSuite)
{
    const auto parameters =
    R"(
    {
        "solution_type" : "Quasi-Static"
    }
    )";

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(SchemeFactoryType::Create(Parameters{parameters}),
                                      "scheme_type is not defined, aborting")
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_Throws_WhenSolutionTypeIsUndefined, KratosGeoMechanicsFastSuite)
{
    const auto parameters =
            R"(
    {
        "scheme_type" : "Backward_Euler"
    }
    )";

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(SchemeFactoryType::Create(Parameters{parameters}),
                                      "solution_type is not defined, aborting")
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_Throws_WhenSchemeTypeOrSolutionTypeAreNotSupported, KratosGeoMechanicsFastSuite)
{
    const auto parameters =
            R"(
    {
        "scheme_type" : "NotSupported",
        "solution_type" : "NotSupported"
    }
    )";

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(SchemeFactoryType::Create(Parameters{parameters}),
                                      "Specified solution_type/scheme_type is not supported, aborting")
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForBackwardEulerQuasiStatic, KratosGeoMechanicsFastSuite)
{
    const auto parameters =
            R"(
    {
        "scheme_type" : "Backward_Euler",
        "solution_type" : "Quasi-Static"
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto backward_euler_scheme = dynamic_cast<const BackwardEulerQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType>*>(scheme.get());

    KRATOS_EXPECT_NE(backward_euler_scheme, nullptr);
}

}

