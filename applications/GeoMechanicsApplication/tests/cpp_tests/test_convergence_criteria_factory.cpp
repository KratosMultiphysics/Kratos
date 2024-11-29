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

#include "custom_utilities/convergence_criteria_factory.hpp"
#include "geo_mechanics_fast_suite.h"
#include "spaces/ublas_space.h"

using namespace Kratos;
using SparseSpaceType                = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType                 = UblasSpace<double, Matrix, Vector>;
using ConvergenceCriteriaFactoryType = ConvergenceCriteriaFactory<SparseSpaceType, LocalSpaceType>;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(Create_ReturnsCorrectConvergenceCriteria_ForDisplacement, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::string valid_parameters = R"(
    {
        "convergence_criterion":              "displacement_criterion",
        "displacement_relative_tolerance":    1.0E-4,
        "displacement_absolute_tolerance":    1.0E-9
    }
    )";

    const auto convergence_criteria = ConvergenceCriteriaFactoryType::Create(Parameters{valid_parameters});
    const auto displacement_criterion =
        dynamic_cast<const DisplacementCriteria<SparseSpaceType, LocalSpaceType>*>(
            convergence_criteria.get());
    KRATOS_EXPECT_NE(displacement_criterion, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(Create_ReturnsCorrectConvergenceCriteria_ForResidual, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::string valid_parameters = R"(
    {
        "convergence_criterion":          "residual_criterion",
        "residual_relative_tolerance":    1.0E-4,
        "residual_absolute_tolerance":    1.0E-9
    }
    )";

    const auto convergence_criteria = ConvergenceCriteriaFactoryType::Create(Parameters{valid_parameters});
    const auto residual_criterion =
        dynamic_cast<const ResidualCriteria<SparseSpaceType, LocalSpaceType>*>(convergence_criteria.get());
    KRATOS_EXPECT_NE(residual_criterion, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(Create_Throws_WhenConvergenceCriterionDoesNotExist,
                          KratosGeoMechanicsFastSuiteWithoutKernel){
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConvergenceCriteriaFactoryType::Create(Parameters{"{}"}),
                                      "No convergence_criterion is defined, aborting.")}

KRATOS_TEST_CASE_IN_SUITE(Create_Throws_WhenConvergenceCriterionIsUnknown, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto invalid_parameters = Parameters{R"({"convergence_criterion" : "something_unknown" })"};
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ConvergenceCriteriaFactoryType::Create(invalid_parameters),
                                      "The convergence_criterion (something_unknown) is unknown, "
                                      "supported criteria are: 'displacement_criterion'")
}

} // namespace Kratos::Testing
