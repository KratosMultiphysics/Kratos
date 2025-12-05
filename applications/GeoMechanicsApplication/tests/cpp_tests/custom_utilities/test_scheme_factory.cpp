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

#include "custom_utilities/scheme_factory.hpp"
#include "includes/expect.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

using namespace Kratos;

namespace Kratos::Testing
{

using SparseSpaceType   = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType    = UblasSpace<double, Matrix, Vector>;
using SchemeFactoryType = SchemeFactory<SparseSpaceType, LocalSpaceType>;

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateScheme_Throws_WhenSchemeTypeIsUndefined)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateScheme_Throws_WhenSolutionTypeIsUndefined)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateScheme_Throws_WhenSchemeTypeOrSolutionTypeAreNotSupported)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "NotSupported",
        "solution_type" : "NotSupported"
    }
    )";

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(SchemeFactoryType::Create(Parameters{parameters}),
                                      "Specified combination of solution_type (NotSupported) and "
                                      "scheme_type (NotSupported) is not supported, aborting")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateScheme_ReturnsCorrectScheme_ForBackwardEulerQuasiStatic)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "Backward_Euler",
        "solution_type" : "Quasi-Static"
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto backward_euler_scheme =
        dynamic_cast<const BackwardEulerQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(backward_euler_scheme, nullptr);
}

} // namespace Kratos::Testing
