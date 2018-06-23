//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "includes/serializer.h"
#include "testing/testing.h"
#include "custom_mappers/nearest_element_mapper.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos {
namespace Testing {

using MappingWeightsVector = std::vector<double>;
using EquationIdVectorType = std::vector<std::size_t>;

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    const Point coords(1.0, 2.45, -23.8);

    const std::size_t source_local_sys_idx = 123;
    const std::size_t dummy_rank = 78;

    NearestElementInterfaceInfo nearest_element_info(coords, source_local_sys_idx, 0);

    const auto nearest_element_info_1(nearest_element_info.Create());
    const auto nearest_element_info_2(nearest_element_info.Create(coords, source_local_sys_idx));
    const auto nearest_element_info_3(nearest_element_info.Create(coords, source_local_sys_idx, dummy_rank));

    // Test if the "Create" function returns the correct object
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(*nearest_element_info_1));
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(*nearest_element_info_2));
    KRATOS_CHECK_EQUAL(typeid(nearest_element_info), typeid(*nearest_element_info_3));
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_ValidProjectionsExist, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_Approximation, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementInterfaceInfo_Serialization, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_BasicTests, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystem, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

KRATOS_TEST_CASE_IN_SUITE(NearestElementLocalSystem_ComputeLocalSystemWithApproximation, KratosMappingApplicationSerialTestSuite)
{
    KRATOS_CHECK(false); // TODO implement test!
}

}  // namespace Testing
}  // namespace Kratos