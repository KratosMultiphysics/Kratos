//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "integration/quadrilateral_gauss_lobatto_integration_points.h"
#include "integration/hexahedron_gauss_lobatto_integration_points.h"
#include "integration/triangle_gauss_lobatto_integration_points.h"
#include "integration/tetrahedron_gauss_lobatto_integration_points.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GaussLobattoQuadrilateralQuadraturesTest, KratosCoreFastSuite)
{

    // In This test we evaluate the Gauss-Lobatto quadratures for integrating
    // f = (x+y) and g = (x+y)^2 over a [-1, 1] quadrilateral

    const auto& r_lobatto_1 = QuadrilateralGaussLobattoIntegrationPoints1();
    const auto& r_lobatto_2 = QuadrilateralGaussLobattoIntegrationPoints2();
    const auto& r_lobatto_3 = QuadrilateralGaussLobattoIntegrationPoints3();

    // Analytical results, reference
    const double integral_f = 0.0;
    const double integral_g = 8.0 / 3.0;

    double quadrature_integral_f = 0.0;
    double quadrature_integral_g = 0.0;

    // Integral for f with Lobatto 1
    for (IndexType IP = 0; IP < r_lobatto_1.IntegrationPoints().size(); ++IP) {
        const auto& r_IP = r_lobatto_1.IntegrationPoints()[IP];
        const double X = r_IP.X();
        const double Y = r_IP.Y();
        const double w = r_IP.Weight();

        quadrature_integral_f += w * (X + Y);
    }

    KRATOS_CHECK_NEAR(quadrature_integral_f, integral_f, 1.0E-6);
    quadrature_integral_f = 0.0;

    // Integral for f with Lobatto 2
    for (IndexType IP = 0; IP < r_lobatto_2.IntegrationPoints().size(); ++IP) {
        const auto& r_IP = r_lobatto_2.IntegrationPoints()[IP];
        const double X = r_IP.X();
        const double Y = r_IP.Y();
        const double w = r_IP.Weight();

        quadrature_integral_f += w * (X + Y);
    }
    KRATOS_CHECK_NEAR(quadrature_integral_f, integral_f, 1.0E-6);

    // Integral for g with Lobatto 3
    for (IndexType IP = 0; IP < r_lobatto_3.IntegrationPoints().size(); ++IP) {
        const auto& r_IP = r_lobatto_3.IntegrationPoints()[IP];
        const double X = r_IP.X();
        const double Y = r_IP.Y();
        const double w = r_IP.Weight();

        quadrature_integral_g += w * (X + Y) * (X + Y);
    }
    KRATOS_CHECK_NEAR(quadrature_integral_g, integral_g, 1.0E-6);

}

KRATOS_TEST_CASE_IN_SUITE(GaussLobattoHexaQuadraturesTest, KratosCoreFastSuite)
{

    // In This test we evaluate the Gauss-Lobatto quadratures for integrating
    // f = (x+y+z+5) over a [-1, 1] hexa

    const auto& r_lobatto = HexahedronGaussLobattoIntegrationPoints2();

    // Analytical results, reference
    const double integral_f = 40.0;

    double quadrature_integral_f = 0.0;

    // Integral for f with Lobatto 1
    for (IndexType IP = 0; IP < r_lobatto.IntegrationPoints().size(); ++IP) {
        const auto& r_IP = r_lobatto.IntegrationPoints()[IP];
        const double X = r_IP.X();
        const double Y = r_IP.Y();
        const double Z = r_IP.Z();
        const double w = r_IP.Weight();

        quadrature_integral_f += w * (X + Y + Z + 5.0);
    }

    KRATOS_CHECK_NEAR(quadrature_integral_f, integral_f, 1.0E-6);
}

KRATOS_TEST_CASE_IN_SUITE(GaussLobattoTriangleQuadraturesTest, KratosCoreFastSuite)
{

    // In This test we evaluate the Gauss-Lobatto quadrature for integrating
    // f = (x+y) over a [0, 1]x[0, 1-x] triangle

    const auto& r_lobatto_1 = TriangleGaussLobattoIntegrationPoints1();

    // Analytical result, reference
    const double integral_f = 1.0 / 3.0;

    double quadrature_integral_f = 0.0;

    // Integral for f with Lobatto 1
    for (IndexType IP = 0; IP < r_lobatto_1.IntegrationPoints().size(); ++IP) {
        const auto& r_IP = r_lobatto_1.IntegrationPoints()[IP];
        const double X = r_IP.X();
        const double Y = r_IP.Y();
        const double w = r_IP.Weight();

        quadrature_integral_f += w * (X + Y);
    }

    KRATOS_CHECK_NEAR(quadrature_integral_f, integral_f, 1.0E-6);
}

KRATOS_TEST_CASE_IN_SUITE(GaussLobattoTetrahedronQuadraturesTest, KratosCoreFastSuite)
{

    // In This test we evaluate the Gauss-Lobatto quadrature for integrating
    // f = 1.0 over an isoparametric tetrahedron

    const auto& r_lobatto_1 = TetrahedronGaussLobattoIntegrationPoints1();

    // Analytical result, reference
    const double integral_f = 1.0 / 6.0;

    double quadrature_integral_f = 0.0;

    // Integral for f with Lobatto 1
    for (IndexType IP = 0; IP < r_lobatto_1.IntegrationPoints().size(); ++IP) {
        const double w = r_lobatto_1.IntegrationPoints()[IP].Weight();
        quadrature_integral_f += w;
    }

    KRATOS_CHECK_NEAR(quadrature_integral_f, integral_f, 1.0E-6);
}

} // namespace Kratos::Testing