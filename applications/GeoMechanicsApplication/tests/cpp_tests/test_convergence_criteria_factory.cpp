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
#include "custom_utilities/convergence_criteria_factory.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"

const std::string validParameters = R"(
    {
        "convergence_criterion":              "displacement_criterion",
        "displacement_relative_tolerance":    1.0E-4,
        "displacement_absolute_tolerance":    1.0E-9
    }
)";


using namespace Kratos;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(Create_ReturnsCorrectConvergenceCriteria_ForDisplacement, KratosGeoMechanicsFastSuite)
{
    const auto convergence_criteria = ConvergenceCriteriaFactory<SparseSpaceType, LocalSpaceType>::Create(Parameters{validParameters});
    const auto displacementCriterion = dynamic_cast<DisplacementCriteria<SparseSpaceType, LocalSpaceType>*>(convergence_criteria.get());
    KRATOS_EXPECT_NE(displacementCriterion, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(Create_Throws_WhenConvergenceCriterionDoesNotExist, KratosGeoMechanicsFastSuite)
{
    ConvergenceCriteriaFactory<SparseSpaceType, LocalSpaceType> factory;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(factory.Create(Parameters{"{}"}), "No convergence_criterion is defined, aborting.")
}

KRATOS_TEST_CASE_IN_SUITE(Create_Throws_WhenConvergenceCriterionIsUnknown, KratosGeoMechanicsFastSuite)
{
    ConvergenceCriteriaFactory<SparseSpaceType, LocalSpaceType> factory;

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(factory.Create(Parameters{R"({"convergence_criterion" : "something_unknown" })"}),
                                      "The convergence_criterion (something_unknown) is unknown, supported criteria are: 'displacement_criterion'")
}



}

