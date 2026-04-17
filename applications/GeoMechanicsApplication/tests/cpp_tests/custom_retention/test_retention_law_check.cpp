// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_retention/retention_law.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace ::testing;

namespace Kratos
{
class MockRetentionLaw : public RetentionLaw
{
public:
    MOCK_METHOD(RetentionLaw::Pointer, Clone, (), (const, override));
    MOCK_METHOD(double, CalculateSaturation, (Parameters & rParameters), (const, override));
    MOCK_METHOD(double, CalculateEffectiveSaturation, (Parameters & rParameters), (const, override));
    MOCK_METHOD(double, CalculateDerivativeOfSaturation, (Parameters & rParameters), (const, override));
    MOCK_METHOD(double, CalculateRelativePermeability, (Parameters & rParameters), (const, override));
    MOCK_METHOD(double, CalculateBishopCoefficient, (Parameters & rParameters), (const, override));
    MOCK_METHOD(void, PrintInfo, (std::ostream & rOStream), (const, override));
    MOCK_METHOD(void, PrintData, (std::ostream & rOStream), (const, override));
    MOCK_METHOD(int, Check, (const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo), (override));
    MOCK_METHOD(void, save, (Serializer & rSerializer), (const, override));
    MOCK_METHOD(void, load, (Serializer & rSerializer), (override));
};
} // namespace Kratos

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(RetentionLaw_CheckEmptyRetentionLawVector, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    std::vector<RetentionLaw::Pointer> retention_law_vector;

    // Act and Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(RetentionLaw::Check(retention_law_vector, Properties{}, ProcessInfo{}),
                                      "A retention law has to be provided.")
}

KRATOS_TEST_CASE_IN_SUITE(RetentionLaw_CheckRetentionLawVectorReturnsOne, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    std::vector<RetentionLaw::Pointer> retention_law_vector;
    const auto                         mock_retention_law = std::make_shared<MockRetentionLaw>();
    retention_law_vector.emplace_back(mock_retention_law);
    retention_law_vector.emplace_back(mock_retention_law);
    EXPECT_CALL(*mock_retention_law, Check(_, _)).WillOnce(Return(1));

    // Act and Assert
    KRATOS_EXPECT_EQ(RetentionLaw::Check(retention_law_vector, Properties{}, ProcessInfo{}), 1);
}

KRATOS_TEST_CASE_IN_SUITE(RetentionLaw_CheckRetentionLawVectorReturnsZero, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    std::vector<RetentionLaw::Pointer> retention_law_vector;
    const auto                         mock_retention_law = std::make_shared<MockRetentionLaw>();
    retention_law_vector.emplace_back(mock_retention_law);
    retention_law_vector.emplace_back(mock_retention_law);
    EXPECT_CALL(*mock_retention_law, Check(_, _)).WillOnce(Return(0));

    // Act and Assert
    KRATOS_EXPECT_EQ(RetentionLaw::Check(retention_law_vector, Properties{}, ProcessInfo{}), 0);
}
} // namespace Kratos::Testing
