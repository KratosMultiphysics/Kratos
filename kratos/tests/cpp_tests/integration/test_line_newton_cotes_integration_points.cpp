//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "integration/line_newton_cotes_integration_points.h"

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(LineNewtonCotesIntegrationPoints1Test, KratosCoreFastSuite)
{
    auto& r_integration_points = LineNewtonCotesIntegrationPoints1::IntegrationPoints();
    KRATOS_CHECK_EQUAL(LineNewtonCotesIntegrationPoints1::IntegrationPointsNumber(), 3);
    KRATOS_CHECK_NEAR(r_integration_points[0].X(), -0.5, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[0].Weight(), 1.333333333333333, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].X(), 0.0, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].Weight(), -0.6666666666666665, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].X(), 0.5, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].Weight(), 1.333333333333333, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(LineNewtonCotesIntegrationPoints2Test, KratosCoreFastSuite)
{
    auto& r_integration_points = LineNewtonCotesIntegrationPoints2::IntegrationPoints();
    KRATOS_CHECK_EQUAL(LineNewtonCotesIntegrationPoints2::IntegrationPointsNumber(), 5);
    KRATOS_CHECK_NEAR(r_integration_points[0].X(), -0.6666666666666666, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[0].Weight(), 1.1, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].X(), -0.3333333333333333, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].Weight(), -1.4, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].X(), 0.0, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].Weight(), 2.6, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].X(), 0.3333333333333333, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].Weight(), -1.4, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].X(), 0.6666666666666666, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].Weight(), 1.1, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(LineNewtonCotesIntegrationPoints3Test, KratosCoreFastSuite)
{
    auto& r_integration_points = LineNewtonCotesIntegrationPoints3::IntegrationPoints();
    KRATOS_CHECK_EQUAL(LineNewtonCotesIntegrationPoints3::IntegrationPointsNumber(), 7);
    KRATOS_CHECK_NEAR(r_integration_points[0].X(), -0.75, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[0].Weight(), 0.9735449735449742, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].X(), -0.5, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].Weight(), -2.019047619047615, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].X(), -0.25, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].Weight(), 4.647619047619042, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].X(), 0.0, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].Weight(), -5.204232804232804, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].X(), 0.25, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].Weight(), 4.647619047619049, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[5].X(), 0.5, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[5].Weight(), -2.019047619047616, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[6].X(), 0.75, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[6].Weight(), 0.9735449735449739, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(LineNewtonCotesIntegrationPoints4Test, KratosCoreFastSuite)
{
    auto& r_integration_points = LineNewtonCotesIntegrationPoints4::IntegrationPoints();
    KRATOS_CHECK_EQUAL(LineNewtonCotesIntegrationPoints4::IntegrationPointsNumber(), 9);
    KRATOS_CHECK_NEAR(r_integration_points[0].X(), -0.8, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[0].Weight(), 0.8917548500881828, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].X(), -0.6, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].Weight(), -2.577160493827184, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].X(), -0.4, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].Weight(), 7.350088183421553, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].X(), -0.2, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].Weight(), -12.14065255731907, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].X(), 0.0, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].Weight(), 14.95194003527322, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[5].X(), 0.2, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[5].Weight(), -12.14065255731914, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[6].X(), 0.4, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[6].Weight(), 7.350088183421514, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[7].X(), 0.6, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[7].Weight(), -2.577160493827156, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[8].X(), 0.8, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[8].Weight(), 0.8917548500881831, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(LineNewtonCotesIntegrationPoints5Test, KratosCoreFastSuite)
{
    auto& r_integration_points = LineNewtonCotesIntegrationPoints5::IntegrationPoints();
    KRATOS_CHECK_EQUAL(LineNewtonCotesIntegrationPoints5::IntegrationPointsNumber(), 11);
    KRATOS_CHECK_NEAR(r_integration_points[0].X(), -0.8333333333333334, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[0].Weight(), 0.8334199134199013, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].X(), -0.6666666666666666, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[1].Weight(), -3.097056277056003, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].X(), -0.5, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[2].Weight(), 10.65437229437228, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].X(), -0.3333333333333333, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[3].Weight(), -23.0561038961028, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].X(), -0.1666666666666667, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[4].Weight(), 37.05246753246604, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[5].X(), 0.0, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[5].Weight(), -42.77419913420147, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[6].X(), 0.1666666666666667, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[6].Weight(), 37.05246753246799, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[7].X(), 0.3333333333333333, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[7].Weight(), -23.05610389610366, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[8].X(), 0.5, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[8].Weight(), 10.65437229437222, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[9].X(), 0.6666666666666666, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[9].Weight(), -3.097056277056263, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[10].X(), 0.8333333333333334, 1e-12);
    KRATOS_CHECK_NEAR(r_integration_points[10].Weight(), 0.8334199134199121, 1e-12);
}

} // namespace Kratos::Testing