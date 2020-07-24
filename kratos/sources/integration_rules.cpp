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

//LINE:

//Gauss-Legendre

namespace Kratos
{
LineGaussLegendreIntegrationPoints1::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints1::msIntegrationPoints = LineGaussLegendreIntegrationPoints1::IntegrationPoints();
LineGaussLegendreIntegrationPoints2::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints2::msIntegrationPoints = LineGaussLegendreIntegrationPoints2::IntegrationPoints();
LineGaussLegendreIntegrationPoints3::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints3::msIntegrationPoints = LineGaussLegendreIntegrationPoints3::IntegrationPoints();
LineGaussLegendreIntegrationPoints4::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints4::msIntegrationPoints = LineGaussLegendreIntegrationPoints4::IntegrationPoints();
LineGaussLegendreIntegrationPoints5::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints5::msIntegrationPoints = LineGaussLegendreIntegrationPoints5::IntegrationPoints();
LineGaussLegendreIntegrationPoints6::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints6::msIntegrationPoints = LineGaussLegendreIntegrationPoints6::IntegrationPoints();
LineGaussLegendreIntegrationPoints7::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints7::msIntegrationPoints = LineGaussLegendreIntegrationPoints7::IntegrationPoints();
LineGaussLegendreIntegrationPoints8::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints8::msIntegrationPoints = LineGaussLegendreIntegrationPoints8::IntegrationPoints();
LineGaussLegendreIntegrationPoints9::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints9::msIntegrationPoints = LineGaussLegendreIntegrationPoints9::IntegrationPoints();
LineGaussLegendreIntegrationPoints10::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints10::msIntegrationPoints = LineGaussLegendreIntegrationPoints10::IntegrationPoints();

//Gauss-Lobatto
LineGaussLobattoIntegrationPoints1::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints1::msIntegrationPoints = LineGaussLobattoIntegrationPoints1::IntegrationPoints();
LineGaussLobattoIntegrationPoints2::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints2::msIntegrationPoints = LineGaussLobattoIntegrationPoints2::IntegrationPoints();
LineGaussLobattoIntegrationPoints3::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints3::msIntegrationPoints = LineGaussLobattoIntegrationPoints3::IntegrationPoints();
LineGaussLobattoIntegrationPoints4::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints4::msIntegrationPoints = LineGaussLobattoIntegrationPoints4::IntegrationPoints();
LineGaussLobattoIntegrationPoints5::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints5::msIntegrationPoints = LineGaussLobattoIntegrationPoints5::IntegrationPoints();
LineGaussLobattoIntegrationPoints6::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints6::msIntegrationPoints = LineGaussLobattoIntegrationPoints6::IntegrationPoints();
LineGaussLobattoIntegrationPoints7::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints7::msIntegrationPoints = LineGaussLobattoIntegrationPoints7::IntegrationPoints();
LineGaussLobattoIntegrationPoints8::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints8::msIntegrationPoints = LineGaussLobattoIntegrationPoints8::IntegrationPoints();
LineGaussLobattoIntegrationPoints9::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints9::msIntegrationPoints = LineGaussLobattoIntegrationPoints9::IntegrationPoints();
LineGaussLobattoIntegrationPoints10::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints10::msIntegrationPoints = LineGaussLobattoIntegrationPoints10::IntegrationPoints();

// Collocation
LineCollocationIntegrationPoints1::IntegrationPointsArrayType LineCollocationIntegrationPoints1::msIntegrationPoints = LineCollocationIntegrationPoints1::IntegrationPoints();
LineCollocationIntegrationPoints2::IntegrationPointsArrayType LineCollocationIntegrationPoints2::msIntegrationPoints = LineCollocationIntegrationPoints2::IntegrationPoints();
LineCollocationIntegrationPoints3::IntegrationPointsArrayType LineCollocationIntegrationPoints3::msIntegrationPoints = LineCollocationIntegrationPoints3::IntegrationPoints();
LineCollocationIntegrationPoints4::IntegrationPointsArrayType LineCollocationIntegrationPoints4::msIntegrationPoints = LineCollocationIntegrationPoints4::IntegrationPoints();
LineCollocationIntegrationPoints5::IntegrationPointsArrayType LineCollocationIntegrationPoints5::msIntegrationPoints = LineCollocationIntegrationPoints5::IntegrationPoints();

//TRIANGLE:

//Gauss-Legendre
TriangleGaussLegendreIntegrationPoints1::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints1::msIntegrationPoints = TriangleGaussLegendreIntegrationPoints1::IntegrationPoints();
TriangleGaussLegendreIntegrationPoints2::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints2::msIntegrationPoints = TriangleGaussLegendreIntegrationPoints2::IntegrationPoints();
TriangleGaussLegendreIntegrationPoints3::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints3::msIntegrationPoints = TriangleGaussLegendreIntegrationPoints3::IntegrationPoints();
TriangleGaussLegendreIntegrationPoints4::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints4::msIntegrationPoints = TriangleGaussLegendreIntegrationPoints4::IntegrationPoints();
TriangleGaussLegendreIntegrationPoints5::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints5::msIntegrationPoints = TriangleGaussLegendreIntegrationPoints5::IntegrationPoints();

// Gauss-Radau
TriangleGaussRadauIntegrationPoints1::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints1::msIntegrationPoints = TriangleGaussRadauIntegrationPoints1::IntegrationPoints();
TriangleGaussRadauIntegrationPoints2::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints2::msIntegrationPoints = TriangleGaussRadauIntegrationPoints2::IntegrationPoints();
TriangleGaussRadauIntegrationPoints3::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints3::msIntegrationPoints = TriangleGaussRadauIntegrationPoints3::IntegrationPoints();
TriangleGaussRadauIntegrationPoints4::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints4::msIntegrationPoints = TriangleGaussRadauIntegrationPoints4::IntegrationPoints();
TriangleGaussRadauIntegrationPoints5::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints5::msIntegrationPoints = TriangleGaussRadauIntegrationPoints5::IntegrationPoints();
TriangleGaussRadauIntegrationPoints6::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints6::msIntegrationPoints = TriangleGaussRadauIntegrationPoints6::IntegrationPoints();

// Collocation
TriangleCollocationIntegrationPoints1::IntegrationPointsArrayType TriangleCollocationIntegrationPoints1::msIntegrationPoints = TriangleCollocationIntegrationPoints1::IntegrationPoints();
TriangleCollocationIntegrationPoints2::IntegrationPointsArrayType TriangleCollocationIntegrationPoints2::msIntegrationPoints = TriangleCollocationIntegrationPoints2::IntegrationPoints();
TriangleCollocationIntegrationPoints3::IntegrationPointsArrayType TriangleCollocationIntegrationPoints3::msIntegrationPoints = TriangleCollocationIntegrationPoints3::IntegrationPoints();
TriangleCollocationIntegrationPoints4::IntegrationPointsArrayType TriangleCollocationIntegrationPoints4::msIntegrationPoints = TriangleCollocationIntegrationPoints4::IntegrationPoints();
TriangleCollocationIntegrationPoints5::IntegrationPointsArrayType TriangleCollocationIntegrationPoints5::msIntegrationPoints = TriangleCollocationIntegrationPoints5::IntegrationPoints();

//QUADRILATERAL:

//Gauss-Legendre
QuadrilateralGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints1::msIntegrationPoints = QuadrilateralGaussLegendreIntegrationPoints1::IntegrationPoints();
QuadrilateralGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints2::msIntegrationPoints = QuadrilateralGaussLegendreIntegrationPoints2::IntegrationPoints();
QuadrilateralGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints3::msIntegrationPoints = QuadrilateralGaussLegendreIntegrationPoints3::IntegrationPoints();
QuadrilateralGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints4::msIntegrationPoints = QuadrilateralGaussLegendreIntegrationPoints4::IntegrationPoints();
QuadrilateralGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints5::msIntegrationPoints = QuadrilateralGaussLegendreIntegrationPoints5::IntegrationPoints();

//Gauss-Lobatto
QuadrilateralGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralGaussLobattoIntegrationPoints1::msIntegrationPoints = QuadrilateralGaussLobattoIntegrationPoints1::IntegrationPoints();
QuadrilateralGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralGaussLobattoIntegrationPoints2::msIntegrationPoints = QuadrilateralGaussLobattoIntegrationPoints2::IntegrationPoints();

// Collocation
QuadrilateralCollocationIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints1::msIntegrationPoints = QuadrilateralCollocationIntegrationPoints1::IntegrationPoints();
QuadrilateralCollocationIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints2::msIntegrationPoints = QuadrilateralCollocationIntegrationPoints2::IntegrationPoints();
QuadrilateralCollocationIntegrationPoints3::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints3::msIntegrationPoints = QuadrilateralCollocationIntegrationPoints3::IntegrationPoints();
QuadrilateralCollocationIntegrationPoints4::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints4::msIntegrationPoints = QuadrilateralCollocationIntegrationPoints4::IntegrationPoints();
QuadrilateralCollocationIntegrationPoints5::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints5::msIntegrationPoints = QuadrilateralCollocationIntegrationPoints5::IntegrationPoints();

//TETRAHEDRON:

//Gauss-Legendre
TetrahedronGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints1::msIntegrationPoints = TetrahedronGaussLegendreIntegrationPoints1::IntegrationPoints();
TetrahedronGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints2::msIntegrationPoints = TetrahedronGaussLegendreIntegrationPoints2::IntegrationPoints();
TetrahedronGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints3::msIntegrationPoints = TetrahedronGaussLegendreIntegrationPoints3::IntegrationPoints();
TetrahedronGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints4::msIntegrationPoints = TetrahedronGaussLegendreIntegrationPoints4::IntegrationPoints();
TetrahedronGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints5::msIntegrationPoints = TetrahedronGaussLegendreIntegrationPoints5::IntegrationPoints();

//PRISM:

//Gauss-Legendre
PrismGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints1::msIntegrationPoints = PrismGaussLegendreIntegrationPoints1::IntegrationPoints();
PrismGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints2::msIntegrationPoints = PrismGaussLegendreIntegrationPoints2::IntegrationPoints();
PrismGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints3::msIntegrationPoints = PrismGaussLegendreIntegrationPoints3::IntegrationPoints();
PrismGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints4::msIntegrationPoints = PrismGaussLegendreIntegrationPoints4::IntegrationPoints();
PrismGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints5::msIntegrationPoints = PrismGaussLegendreIntegrationPoints5::IntegrationPoints();

/* EXTENDED VALUES (just one point along the plane) */
PrismGaussLegendreIntegrationPointsExt1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt1::msIntegrationPoints = PrismGaussLegendreIntegrationPointsExt1::IntegrationPoints();
PrismGaussLegendreIntegrationPointsExt2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt2::msIntegrationPoints = PrismGaussLegendreIntegrationPointsExt2::IntegrationPoints();
PrismGaussLegendreIntegrationPointsExt3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt3::msIntegrationPoints = PrismGaussLegendreIntegrationPointsExt3::IntegrationPoints();
PrismGaussLegendreIntegrationPointsExt4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt4::msIntegrationPoints = PrismGaussLegendreIntegrationPointsExt4::IntegrationPoints();
PrismGaussLegendreIntegrationPointsExt5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt5::msIntegrationPoints = PrismGaussLegendreIntegrationPointsExt5::IntegrationPoints();

//Gauss-Lobatto
PrismGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
PrismGaussLobattoIntegrationPoints1::msIntegrationPoints = PrismGaussLobattoIntegrationPoints1::IntegrationPoints();
PrismGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
PrismGaussLobattoIntegrationPoints2::msIntegrationPoints = PrismGaussLobattoIntegrationPoints2::IntegrationPoints();

//HEXAHEDRON:

//Gauss-Legendre
HexahedronGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints1::msIntegrationPoints = HexahedronGaussLegendreIntegrationPoints1::IntegrationPoints();
HexahedronGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints2::msIntegrationPoints = HexahedronGaussLegendreIntegrationPoints2::IntegrationPoints();
HexahedronGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints3::msIntegrationPoints = HexahedronGaussLegendreIntegrationPoints3::IntegrationPoints();
HexahedronGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints4::msIntegrationPoints = HexahedronGaussLegendreIntegrationPoints4::IntegrationPoints();
HexahedronGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints5::msIntegrationPoints = HexahedronGaussLegendreIntegrationPoints5::IntegrationPoints();

//Gauss-Lobatto
HexahedronGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
HexahedronGaussLobattoIntegrationPoints1::msIntegrationPoints = HexahedronGaussLobattoIntegrationPoints1::IntegrationPoints();
HexahedronGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
HexahedronGaussLobattoIntegrationPoints2::msIntegrationPoints = HexahedronGaussLobattoIntegrationPoints2::IntegrationPoints();
}
