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

#include <gtest/gtest.h>
#include <gmock/gmock.h>

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

private:
    MOCK_METHOD(void, save, (Serializer & rSerializer));
    MOCK_METHOD(void, load, (Serializer & rSerializer));
};
}
namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(RetentionLaw_CheckRetentionwVector,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    std::vector<RetentionLaw::Pointer> retention_law_vector;
    const auto                         properties = Properties();
    const auto                         current_process_info = ProcessInfo();

    // Act and Assert
    KRATOS_EXPECT_EQ(RetentionLaw::Check(retention_law_vector, properties, current_process_info), 0);

    // Arrange
    const auto mock_retention_law = std::make_shared<MockRetentionLaw>();
    retention_law_vector.push_back(mock_retention_law);
    retention_law_vector.push_back(mock_retention_law);
    EXPECT_CALL(*mock_retention_law,Check(_,_)).WillOnce(Return(1));

    // Act and Assert
    KRATOS_EXPECT_EQ(RetentionLaw::Check(retention_law_vector, properties, current_process_info), 1);

    EXPECT_CALL(*mock_retention_law, Check(_, _)).WillOnce(Return(0));

    // Act and Assert
    KRATOS_EXPECT_EQ(RetentionLaw::Check(retention_law_vector, properties, current_process_info), 0);
}

} // namespace Kratos::Testing
