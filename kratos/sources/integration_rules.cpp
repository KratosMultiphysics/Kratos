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

//Gauss-Legendre
LineGaussLegendreIntegrationPoints1::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints1::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints1::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints2::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints2::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints2::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints3::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints3::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints3::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints4::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints4::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints4::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints5::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints5::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints5::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints6::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints6::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints6::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints7::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints7::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints7::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints8::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints8::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints8::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints9::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints9::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints9::IntegrationPoints();}();
LineGaussLegendreIntegrationPoints10::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints10::msIntegrationPoints = [](){return LineGaussLegendreIntegrationPoints10::IntegrationPoints();}();

//Gauss-Lobatto
LineGaussLobattoIntegrationPoints1::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints1::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints1::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints2::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints2::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints2::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints3::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints3::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints3::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints4::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints4::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints4::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints5::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints5::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints5::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints6::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints6::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints6::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints7::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints7::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints7::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints8::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints8::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints8::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints9::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints9::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints9::IntegrationPoints();}();
LineGaussLobattoIntegrationPoints10::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints10::msIntegrationPoints = [](){return LineGaussLobattoIntegrationPoints10::IntegrationPoints();}();

// Collocation
LineCollocationIntegrationPoints1::IntegrationPointsArrayType LineCollocationIntegrationPoints1::msIntegrationPoints = [](){return LineCollocationIntegrationPoints1::IntegrationPoints();}();
LineCollocationIntegrationPoints2::IntegrationPointsArrayType LineCollocationIntegrationPoints2::msIntegrationPoints = [](){return LineCollocationIntegrationPoints2::IntegrationPoints();}();
LineCollocationIntegrationPoints3::IntegrationPointsArrayType LineCollocationIntegrationPoints3::msIntegrationPoints = [](){return LineCollocationIntegrationPoints3::IntegrationPoints();}();
LineCollocationIntegrationPoints4::IntegrationPointsArrayType LineCollocationIntegrationPoints4::msIntegrationPoints = [](){return LineCollocationIntegrationPoints4::IntegrationPoints();}();
LineCollocationIntegrationPoints5::IntegrationPointsArrayType LineCollocationIntegrationPoints5::msIntegrationPoints = [](){return LineCollocationIntegrationPoints5::IntegrationPoints();}();

//TRIANGLE:

//Gauss-Legendre
TriangleGaussLegendreIntegrationPoints1::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints1::msIntegrationPoints = [](){return TriangleGaussLegendreIntegrationPoints1::IntegrationPoints();}();
TriangleGaussLegendreIntegrationPoints2::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints2::msIntegrationPoints = [](){return TriangleGaussLegendreIntegrationPoints2::IntegrationPoints();}();
TriangleGaussLegendreIntegrationPoints3::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints3::msIntegrationPoints = [](){return TriangleGaussLegendreIntegrationPoints3::IntegrationPoints();}();
TriangleGaussLegendreIntegrationPoints4::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints4::msIntegrationPoints = [](){return TriangleGaussLegendreIntegrationPoints4::IntegrationPoints();}();
TriangleGaussLegendreIntegrationPoints5::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints5::msIntegrationPoints = [](){return TriangleGaussLegendreIntegrationPoints5::IntegrationPoints();}();

// Gauss-Radau
TriangleGaussRadauIntegrationPoints1::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints1::msIntegrationPoints = [](){return TriangleGaussRadauIntegrationPoints1::IntegrationPoints();}();
TriangleGaussRadauIntegrationPoints2::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints2::msIntegrationPoints = [](){return TriangleGaussRadauIntegrationPoints2::IntegrationPoints();}();
TriangleGaussRadauIntegrationPoints3::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints3::msIntegrationPoints = [](){return TriangleGaussRadauIntegrationPoints3::IntegrationPoints();}();
TriangleGaussRadauIntegrationPoints4::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints4::msIntegrationPoints = [](){return TriangleGaussRadauIntegrationPoints4::IntegrationPoints();}();
TriangleGaussRadauIntegrationPoints5::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints5::msIntegrationPoints = [](){return TriangleGaussRadauIntegrationPoints5::IntegrationPoints();}();
TriangleGaussRadauIntegrationPoints6::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints6::msIntegrationPoints = [](){return TriangleGaussRadauIntegrationPoints6::IntegrationPoints();}();

// Collocation
TriangleCollocationIntegrationPoints1::IntegrationPointsArrayType TriangleCollocationIntegrationPoints1::msIntegrationPoints = [](){return TriangleCollocationIntegrationPoints1::IntegrationPoints();}();
TriangleCollocationIntegrationPoints2::IntegrationPointsArrayType TriangleCollocationIntegrationPoints2::msIntegrationPoints = [](){return TriangleCollocationIntegrationPoints2::IntegrationPoints();}();
TriangleCollocationIntegrationPoints3::IntegrationPointsArrayType TriangleCollocationIntegrationPoints3::msIntegrationPoints = [](){return TriangleCollocationIntegrationPoints3::IntegrationPoints();}();
TriangleCollocationIntegrationPoints4::IntegrationPointsArrayType TriangleCollocationIntegrationPoints4::msIntegrationPoints = [](){return TriangleCollocationIntegrationPoints4::IntegrationPoints();}();
TriangleCollocationIntegrationPoints5::IntegrationPointsArrayType TriangleCollocationIntegrationPoints5::msIntegrationPoints = [](){return TriangleCollocationIntegrationPoints5::IntegrationPoints();}();

//QUADRILATERAL:

//Gauss-Legendre
QuadrilateralGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints1::msIntegrationPoints = [](){return QuadrilateralGaussLegendreIntegrationPoints1::IntegrationPoints();}();
QuadrilateralGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints2::msIntegrationPoints = [](){return QuadrilateralGaussLegendreIntegrationPoints2::IntegrationPoints();}();
QuadrilateralGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints3::msIntegrationPoints = [](){return QuadrilateralGaussLegendreIntegrationPoints3::IntegrationPoints();}();
QuadrilateralGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints4::msIntegrationPoints = [](){return QuadrilateralGaussLegendreIntegrationPoints4::IntegrationPoints();}();
QuadrilateralGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints5::msIntegrationPoints = [](){return QuadrilateralGaussLegendreIntegrationPoints5::IntegrationPoints();}();

//Gauss-Lobatto
QuadrilateralGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralGaussLobattoIntegrationPoints1::msIntegrationPoints = [](){return QuadrilateralGaussLobattoIntegrationPoints1::IntegrationPoints();}();
QuadrilateralGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralGaussLobattoIntegrationPoints2::msIntegrationPoints = [](){return QuadrilateralGaussLobattoIntegrationPoints2::IntegrationPoints();}();

// Collocation
QuadrilateralCollocationIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints1::msIntegrationPoints = [](){return QuadrilateralCollocationIntegrationPoints1::IntegrationPoints();}();
QuadrilateralCollocationIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints2::msIntegrationPoints = [](){return QuadrilateralCollocationIntegrationPoints2::IntegrationPoints();}();
QuadrilateralCollocationIntegrationPoints3::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints3::msIntegrationPoints = [](){return QuadrilateralCollocationIntegrationPoints3::IntegrationPoints();}();
QuadrilateralCollocationIntegrationPoints4::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints4::msIntegrationPoints = [](){return QuadrilateralCollocationIntegrationPoints4::IntegrationPoints();}();
QuadrilateralCollocationIntegrationPoints5::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints5::msIntegrationPoints = [](){return QuadrilateralCollocationIntegrationPoints5::IntegrationPoints();}();

//TETRAHEDRON:

//Gauss-Legendre
TetrahedronGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints1::msIntegrationPoints = [](){return TetrahedronGaussLegendreIntegrationPoints1::IntegrationPoints();}();
TetrahedronGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints2::msIntegrationPoints = [](){return TetrahedronGaussLegendreIntegrationPoints2::IntegrationPoints();}();
TetrahedronGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints3::msIntegrationPoints = [](){return TetrahedronGaussLegendreIntegrationPoints3::IntegrationPoints();}();
TetrahedronGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints4::msIntegrationPoints = [](){return TetrahedronGaussLegendreIntegrationPoints4::IntegrationPoints();}();
TetrahedronGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints5::msIntegrationPoints = [](){return TetrahedronGaussLegendreIntegrationPoints5::IntegrationPoints();}();

//PRISM:

//Gauss-Legendre
PrismGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints1::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPoints1::IntegrationPoints();}();
PrismGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints2::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPoints2::IntegrationPoints();}();
PrismGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints3::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPoints3::IntegrationPoints();}();
PrismGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints4::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPoints4::IntegrationPoints();}();
PrismGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints5::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPoints5::IntegrationPoints();}();

/* EXTENDED VALUES (just one point along the plane) */
PrismGaussLegendreIntegrationPointsExt1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt1::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPointsExt1::IntegrationPoints();}();
PrismGaussLegendreIntegrationPointsExt2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt2::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPointsExt2::IntegrationPoints();}();
PrismGaussLegendreIntegrationPointsExt3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt3::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPointsExt3::IntegrationPoints();}();
PrismGaussLegendreIntegrationPointsExt4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt4::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPointsExt4::IntegrationPoints();}();
PrismGaussLegendreIntegrationPointsExt5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt5::msIntegrationPoints = [](){return PrismGaussLegendreIntegrationPointsExt5::IntegrationPoints();}();

//Gauss-Lobatto
PrismGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
PrismGaussLobattoIntegrationPoints1::msIntegrationPoints = [](){return PrismGaussLobattoIntegrationPoints1::IntegrationPoints();}();
PrismGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
PrismGaussLobattoIntegrationPoints2::msIntegrationPoints = [](){return PrismGaussLobattoIntegrationPoints2::IntegrationPoints();}();

//HEXAHEDRON:

//Gauss-Legendre
HexahedronGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints1::msIntegrationPoints = [](){return HexahedronGaussLegendreIntegrationPoints1::IntegrationPoints();}();
HexahedronGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints2::msIntegrationPoints = [](){return HexahedronGaussLegendreIntegrationPoints2::IntegrationPoints();}();
HexahedronGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints3::msIntegrationPoints = [](){return HexahedronGaussLegendreIntegrationPoints3::IntegrationPoints();}();
HexahedronGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints4::msIntegrationPoints = [](){return HexahedronGaussLegendreIntegrationPoints4::IntegrationPoints();}();
HexahedronGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints5::msIntegrationPoints = [](){return HexahedronGaussLegendreIntegrationPoints5::IntegrationPoints();}();

//Gauss-Lobatto
HexahedronGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
HexahedronGaussLobattoIntegrationPoints1::msIntegrationPoints = [](){return HexahedronGaussLobattoIntegrationPoints1::IntegrationPoints();}();
HexahedronGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
HexahedronGaussLobattoIntegrationPoints2::msIntegrationPoints = [](){return HexahedronGaussLobattoIntegrationPoints2::IntegrationPoints();}();
}
