//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Project includes
#include "testing/testing.h"
#include "integration/integration_info.h"
#include "geometries/geometry_data.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(CheckIntegrationInfo, KratosCoreFastSuite) {
    IntegrationInfo integration_info(2, GeometryData::IntegrationMethod::GI_GAUSS_2);

    // Definition and setting
    KRATOS_EXPECT_EQ(integration_info.GetNumberOfIntegrationPointsPerSpan(0), 2);
    KRATOS_EXPECT_EQ(integration_info.GetNumberOfIntegrationPointsPerSpan(1), 2);
    KRATOS_EXPECT_EQ(integration_info.GetQuadratureMethod(0), IntegrationInfo::QuadratureMethod::GAUSS);

    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
    KRATOS_EXPECT_NE(integration_info.GetQuadratureMethod(0), IntegrationInfo::QuadratureMethod::GAUSS);
    KRATOS_EXPECT_EQ(integration_info.GetQuadratureMethod(0), IntegrationInfo::QuadratureMethod::GRID);
}

KRATOS_TEST_CASE_IN_SUITE(CheckIntegrationInfoFlags, KratosCoreFastSuite) {
    IntegrationInfo integration_info(2, GeometryData::IntegrationMethod::GI_GAUSS_2);

    KRATOS_EXPECT_FALSE(integration_info.IsDefined(IntegrationInfo::DO_NOT_CREATE_TESSELLATION_ON_SLAVE));
    integration_info.Set(IntegrationInfo::DO_NOT_CREATE_TESSELLATION_ON_SLAVE, false);
    KRATOS_EXPECT_TRUE(integration_info.IsDefined(IntegrationInfo::DO_NOT_CREATE_TESSELLATION_ON_SLAVE));
    KRATOS_EXPECT_FALSE(integration_info.Is(IntegrationInfo::DO_NOT_CREATE_TESSELLATION_ON_SLAVE));
}

}
}