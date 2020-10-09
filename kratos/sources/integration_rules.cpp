//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "integration/line_gauss_legendre_integration_points.h"
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "integration/triangle_gauss_radau_integration_points.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"
#include "integration/tetrahedron_gauss_legendre_integration_points.h"
#include "integration/hexahedron_gauss_legendre_integration_points.h"
#include "integration/prism_gauss_legendre_integration_points.h"

#include "integration/line_gauss_lobatto_integration_points.h"
#include "integration/quadrilateral_gauss_lobatto_integration_points.h"
#include "integration/hexahedron_gauss_lobatto_integration_points.h"
#include "integration/prism_gauss_lobatto_integration_points.h"

#include "integration/line_collocation_integration_points.h"
#include "integration/triangle_collocation_integration_points.h"
#include "integration/quadrilateral_collocation_integration_points.h"

namespace Kratos
{
//LINE:

// Gauss-Legendre
const auto s_integration_points_LineGaussLegendreIntegrationPoints1 = LineGaussLegendreIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints2 = LineGaussLegendreIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints3 = LineGaussLegendreIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints4 = LineGaussLegendreIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints5 = LineGaussLegendreIntegrationPoints5::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints6 = LineGaussLegendreIntegrationPoints6::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints7 = LineGaussLegendreIntegrationPoints7::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints8 = LineGaussLegendreIntegrationPoints8::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints9 = LineGaussLegendreIntegrationPoints9::IntegrationPoints();
const auto s_integration_points_LineGaussLegendreIntegrationPoints10 = LineGaussLegendreIntegrationPoints10::IntegrationPoints();

// Gauss-Lobatto
const auto s_integration_points_LineGaussLobattoIntegrationPoints1 = LineGaussLobattoIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints2 = LineGaussLobattoIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints3 = LineGaussLobattoIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints4 = LineGaussLobattoIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints5 = LineGaussLobattoIntegrationPoints5::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints6 = LineGaussLobattoIntegrationPoints6::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints7 = LineGaussLobattoIntegrationPoints7::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints8 = LineGaussLobattoIntegrationPoints8::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints9 = LineGaussLobattoIntegrationPoints9::IntegrationPoints();
const auto s_integration_points_LineGaussLobattoIntegrationPoints10 = LineGaussLobattoIntegrationPoints10::IntegrationPoints();

// Collocation
const auto s_integration_points_LineCollocationIntegrationPoints1 = LineCollocationIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_LineCollocationIntegrationPoints2 = LineCollocationIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_LineCollocationIntegrationPoints3 = LineCollocationIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_LineCollocationIntegrationPoints4 = LineCollocationIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_LineCollocationIntegrationPoints5 = LineCollocationIntegrationPoints5::IntegrationPoints();

// TRIANGLE:

// Gauss-Legendre
const auto s_integration_points_TriangleGaussLegendreIntegrationPoints1 = TriangleGaussLegendreIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_TriangleGaussLegendreIntegrationPoints2 = TriangleGaussLegendreIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_TriangleGaussLegendreIntegrationPoints3 = TriangleGaussLegendreIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_TriangleGaussLegendreIntegrationPoints4 = TriangleGaussLegendreIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_TriangleGaussLegendreIntegrationPoints5 = TriangleGaussLegendreIntegrationPoints5::IntegrationPoints();

// Gauss-Radau
const auto s_integration_points_TriangleGaussRadauIntegrationPoints1 = TriangleGaussRadauIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_TriangleGaussRadauIntegrationPoints2 = TriangleGaussRadauIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_TriangleGaussRadauIntegrationPoints3 = TriangleGaussRadauIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_TriangleGaussRadauIntegrationPoints4 = TriangleGaussRadauIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_TriangleGaussRadauIntegrationPoints5 = TriangleGaussRadauIntegrationPoints5::IntegrationPoints();
const auto s_integration_points_TriangleGaussRadauIntegrationPoints6 = TriangleGaussRadauIntegrationPoints6::IntegrationPoints();

// Collocation
const auto s_integration_points_TriangleCollocationIntegrationPoints1 = TriangleCollocationIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_TriangleCollocationIntegrationPoints2 = TriangleCollocationIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_TriangleCollocationIntegrationPoints3 = TriangleCollocationIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_TriangleCollocationIntegrationPoints4 = TriangleCollocationIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_TriangleCollocationIntegrationPoints5 = TriangleCollocationIntegrationPoints5::IntegrationPoints();

// QUADRILATERAL:

// Gauss-Legendre
const auto s_integration_points_QuadrilateralGaussLegendreIntegrationPoints1 = QuadrilateralGaussLegendreIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_QuadrilateralGaussLegendreIntegrationPoints2 = QuadrilateralGaussLegendreIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_QuadrilateralGaussLegendreIntegrationPoints3 = QuadrilateralGaussLegendreIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_QuadrilateralGaussLegendreIntegrationPoints4 = QuadrilateralGaussLegendreIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_QuadrilateralGaussLegendreIntegrationPoints5 = QuadrilateralGaussLegendreIntegrationPoints5::IntegrationPoints();

// Gauss-Lobatto
const auto s_integration_points_QuadrilateralGaussLobattoIntegrationPoints1 = QuadrilateralGaussLobattoIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_QuadrilateralGaussLobattoIntegrationPoints2 = QuadrilateralGaussLobattoIntegrationPoints2::IntegrationPoints();

// Collocation
const auto s_integration_points_QuadrilateralCollocationIntegrationPoints1 = QuadrilateralCollocationIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_QuadrilateralCollocationIntegrationPoints2 = QuadrilateralCollocationIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_QuadrilateralCollocationIntegrationPoints3 = QuadrilateralCollocationIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_QuadrilateralCollocationIntegrationPoints4 = QuadrilateralCollocationIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_QuadrilateralCollocationIntegrationPoints5 = QuadrilateralCollocationIntegrationPoints5::IntegrationPoints();

// TETRAHEDRON:

// Gauss-Legendre
const auto s_integration_points_TetrahedronGaussLegendreIntegrationPoints1 = TetrahedronGaussLegendreIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_TetrahedronGaussLegendreIntegrationPoints2 = TetrahedronGaussLegendreIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_TetrahedronGaussLegendreIntegrationPoints3 = TetrahedronGaussLegendreIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_TetrahedronGaussLegendreIntegrationPoints4 = TetrahedronGaussLegendreIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_TetrahedronGaussLegendreIntegrationPoints5 = TetrahedronGaussLegendreIntegrationPoints5::IntegrationPoints();

// PRISM:

// Gauss-Legendre
const auto s_integration_points_PrismGaussLegendreIntegrationPoints1 = PrismGaussLegendreIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPoints2 = PrismGaussLegendreIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPoints3 = PrismGaussLegendreIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPoints4 = PrismGaussLegendreIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPoints5 = PrismGaussLegendreIntegrationPoints5::IntegrationPoints();

/* EXTENDED VALUES (just one point along the plane) */
const auto s_integration_points_PrismGaussLegendreIntegrationPointsExt1 = PrismGaussLegendreIntegrationPointsExt1::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPointsExt2 = PrismGaussLegendreIntegrationPointsExt2::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPointsExt3 = PrismGaussLegendreIntegrationPointsExt3::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPointsExt4 = PrismGaussLegendreIntegrationPointsExt4::IntegrationPoints();
const auto s_integration_points_PrismGaussLegendreIntegrationPointsExt5 = PrismGaussLegendreIntegrationPointsExt5::IntegrationPoints();

// Gauss-Lobatto
const auto s_integration_points_PrismGaussLobattoIntegrationPoints1 = PrismGaussLobattoIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_PrismGaussLobattoIntegrationPoints2 = PrismGaussLobattoIntegrationPoints2::IntegrationPoints();

// HEXAHEDRON:

// Gauss-Legendre
const auto s_integration_points_HexahedronGaussLegendreIntegrationPoints1 = HexahedronGaussLegendreIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_HexahedronGaussLegendreIntegrationPoints2 = HexahedronGaussLegendreIntegrationPoints2::IntegrationPoints();
const auto s_integration_points_HexahedronGaussLegendreIntegrationPoints3 = HexahedronGaussLegendreIntegrationPoints3::IntegrationPoints();
const auto s_integration_points_HexahedronGaussLegendreIntegrationPoints4 = HexahedronGaussLegendreIntegrationPoints4::IntegrationPoints();
const auto s_integration_points_HexahedronGaussLegendreIntegrationPoints5 = HexahedronGaussLegendreIntegrationPoints5::IntegrationPoints();

// Gauss-Lobatto
const auto s_integration_points_HexahedronGaussLobattoIntegrationPoints1 = HexahedronGaussLobattoIntegrationPoints1::IntegrationPoints();
const auto s_integration_points_HexahedronGaussLobattoIntegrationPoints2 = HexahedronGaussLobattoIntegrationPoints2::IntegrationPoints();
}
