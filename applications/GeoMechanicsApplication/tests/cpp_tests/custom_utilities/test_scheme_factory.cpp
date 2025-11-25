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
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

using SparseSpaceType   = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType    = UblasSpace<double, Matrix, Vector>;
using SchemeFactoryType = SchemeFactory<SparseSpaceType, LocalSpaceType>;

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_Throws_WhenSchemeTypeIsUndefined, KratosGeoMechanicsFastSuiteWithoutKernel)
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

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_Throws_WhenSolutionTypeIsUndefined, KratosGeoMechanicsFastSuiteWithoutKernel)
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

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_Throws_WhenSchemeTypeOrSolutionTypeAreNotSupported,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
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

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForBackwardEulerQuasiStatic,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
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

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForBackwardEulerQuasiStaticPw,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "Backward_Euler",
        "solution_type" : "Quasi-Static",
        "solver_type" : "Pw"
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto backward_euler_scheme =
        dynamic_cast<const BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(backward_euler_scheme, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForNewmarkDynamic,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "Newmark",
        "solution_type" : "dynamic",
        "newmark_beta" : 0.25,
        "newmark_gamma" : 0.5,
        "newmark_theta" : 0.5
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto newmark_scheme =
        dynamic_cast<const NewmarkDynamicUPwScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(newmark_scheme, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForNewmarkQuasiStaticUndamped,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "Newmark",
        "solution_type" : "Quasi-Static",
        "newmark_beta" : 0.25,
        "newmark_gamma" : 0.5,
        "newmark_theta" : 0.5
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto newmark_scheme =
        dynamic_cast<const NewmarkQuasistaticUPwScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(newmark_scheme, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForNewmarkQuasiStaticDamped,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "Newmark",
        "solution_type" : "Quasi-Static",
        "newmark_beta" : 0.25,
        "newmark_gamma" : 0.5,
        "newmark_theta" : 0.5,
        "rayleigh_m" : 0.1
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto newmark_scheme =
        dynamic_cast<const NewmarkQuasistaticDampedUPwScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(newmark_scheme, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForNewmarkQuasiStaticDampedWithK,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "Newmark",
        "solution_type" : "Quasi-Static",
        "newmark_beta" : 0.25,
        "newmark_gamma" : 0.5,
        "newmark_theta" : 0.5,
        "rayleigh_k" : 1.0e-5
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto newmark_scheme =
        dynamic_cast<const NewmarkQuasistaticDampedUPwScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(newmark_scheme, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForNewmarkQuasiStaticPw,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "Newmark",
        "solution_type" : "Quasi-Static",
        "newmark_beta" : 0.25,
        "newmark_gamma" : 0.5,
        "newmark_theta" : 0.5,
        "solver_type" : "Pw"
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto newmark_scheme =
        dynamic_cast<const NewmarkQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(newmark_scheme, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForStaticLoadStepping,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "load_stepping",
        "solution_type" : "static"
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto load_stepping_scheme =
        dynamic_cast<const LoadSteppingScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(load_stepping_scheme, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateScheme_ReturnsCorrectScheme_ForStaticDefault,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto parameters =
        R"(
    {
        "scheme_type" : "some_scheme",
        "solution_type" : "static"
    }
    )";

    const auto scheme = SchemeFactoryType::Create(Parameters{parameters});
    const auto static_scheme =
        dynamic_cast<const GeoMechanicsStaticScheme<SparseSpaceType, LocalSpaceType>*>(
            scheme.get());

    KRATOS_EXPECT_NE(static_scheme, nullptr);
}

} // namespace Kratos::Testing
