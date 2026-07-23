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
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <string>

using namespace Kratos;
using SparseSpaceType                = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType                 = UblasSpace<double, Matrix, Vector>;
using ConvergenceCriteriaFactoryType = ConvergenceCriteriaFactory<SparseSpaceType, LocalSpaceType>;
using DisplacementCriterionType      = DisplacementCriteria<SparseSpaceType, LocalSpaceType>;
using ResidualCriterionType          = ResidualCriteria<SparseSpaceType, LocalSpaceType>;
using AndCriterionType               = And_Criteria<SparseSpaceType, LocalSpaceType>;
using OrCriterionType                = Or_Criteria<SparseSpaceType, LocalSpaceType>;
using WaterPressureCriterionType     = MixedGenericCriteria<SparseSpaceType, LocalSpaceType>;
using DisplacementAndWaterPressureCriterionType = AndCriterionType;

using namespace std::string_literals;

namespace Kratos::Testing
{

struct CreateDisplacementCriterionTest {
    using CriterionType = DisplacementCriterionType;
    static const std::string CriterionDefinition;
};

const auto CreateDisplacementCriterionTest::CriterionDefinition = R"(
{
    "convergence_criterion": "displacement_criterion",
    "displacement_relative_tolerance": 1.0E-4,
    "displacement_absolute_tolerance": 1.0E-9
}
)"s;

struct CreateResidualCriterionTest {
    using CriterionType = ResidualCriterionType;
    static const std::string CriterionDefinition;
};

const auto CreateResidualCriterionTest::CriterionDefinition = R"(
{
    "convergence_criterion": "residual_criterion",
    "residual_relative_tolerance": 1.0E-4,
    "residual_absolute_tolerance": 1.0E-9
}
)"s;

struct CreateAndCriterionTest {
    using CriterionType = AndCriterionType;
    static const std::string CriterionDefinition;
};

const auto CreateAndCriterionTest::CriterionDefinition = R"(
{
    "convergence_criterion": "and_criterion",
    "displacement_relative_tolerance": 1.0E-4,
    "displacement_absolute_tolerance": 1.0E-9,
    "residual_relative_tolerance": 1.0E-4,
    "residual_absolute_tolerance": 1.0E-9
}
)"s;

struct CreateOrCriterionTest {
    using CriterionType = OrCriterionType;
    static const std::string CriterionDefinition;
};

const auto CreateOrCriterionTest::CriterionDefinition = R"(
{
    "convergence_criterion": "or_criterion",
    "displacement_relative_tolerance": 1.0E-4,
    "displacement_absolute_tolerance": 1.0E-9,
    "residual_relative_tolerance": 1.0E-4,
    "residual_absolute_tolerance": 1.0E-9
}
)"s;

struct CreateWaterPressureCriterionTest {
    using CriterionType = WaterPressureCriterionType;
    static const std::string CriterionDefinition;
};

const auto CreateWaterPressureCriterionTest::CriterionDefinition = R"(
{
    "convergence_criterion": "water_pressure_criterion",
    "water_pressure_relative_tolerance": 1.0E-4,
    "water_pressure_absolute_tolerance": 1.0E-9
}
)"s;

struct CreateDisplacementAndWaterPressureCriterionTest {
    using CriterionType = DisplacementAndWaterPressureCriterionType;
    static const std::string CriterionDefinition;
};

const auto CreateDisplacementAndWaterPressureCriterionTest::CriterionDefinition = R"(
{
    "convergence_criterion": "displacement_and_water_pressure_criterion",
    "displacement_relative_tolerance": 1.0E-4,
    "displacement_absolute_tolerance": 1.0E-9,
    "water_pressure_relative_tolerance": 1.0E-4,
    "water_pressure_absolute_tolerance": 1.0E-9
}
)"s;

template <typename T>
class CreateConvergenceCriterionTest : public testing::Test
{
};

using CreateConvergenceCriterionTestTypes =
    ::testing::Types<CreateDisplacementCriterionTest, CreateResidualCriterionTest, CreateAndCriterionTest, CreateOrCriterionTest, CreateWaterPressureCriterionTest, CreateDisplacementAndWaterPressureCriterionTest>;
TYPED_TEST_SUITE(CreateConvergenceCriterionTest, CreateConvergenceCriterionTestTypes);

TYPED_TEST(CreateConvergenceCriterionTest, ConvergenceCriteriaFactoryProducesDefinedCriterionType)
{
    const auto p_convergence_criterion =
        ConvergenceCriteriaFactoryType::Create(Parameters{TypeParam::CriterionDefinition});
    KRATOS_EXPECT_NE(dynamic_cast<const TypeParam::CriterionType*>(p_convergence_criterion.get()), nullptr);
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
