// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "integration/line_gauss_legendre_integration_points.h"
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "integration/quadrilateral_gauss_legendre_integration_points.h"
#include "integration/tetrahedron_gauss_legendre_integration_points.h"
#include "integration/hexahedron_gauss_legendre_integration_points.h"
#include "integration/prism_gauss_legendre_integration_points.h"

#include "integration/line_gauss_lobatto_integration_points.h"
#include "integration/quadrilateral_gauss_lobatto_integration_points.h"
#include "integration/hexahedron_gauss_lobatto_integration_points.h"
#include "integration/prism_gauss_lobatto_integration_points.h"

#define tet10_a 0.108103018168070
#define tet10_b 0.445948490915965
#define tet10_c 0.816847572980459

//LINE:

//Gauss-Legendre

namespace Kratos
{
LineGaussLegendreIntegrationPoints1::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType(0.00, 2.00)
    }
};

LineGaussLegendreIntegrationPoints2::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType(-std::sqrt(1.00 / 3.00), 1.00),
        IntegrationPointType( std::sqrt(1.00 / 3.00), 1.00)
    }
};

LineGaussLegendreIntegrationPoints3::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType(-std::sqrt(3.00 / 5.00), 5.00 / 9.00),
        IntegrationPointType( 0.00                  , 8.00 / 9.00),
        IntegrationPointType( std::sqrt(3.00 / 5.00), 5.00 / 9.00)
    }
};

LineGaussLegendreIntegrationPoints4::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.861136311594053, 0.347854845137454),
        IntegrationPointType(-0.339981043584856, 0.652145154862546),
        IntegrationPointType( 0.339981043584856, 0.652145154862546),
        IntegrationPointType( 0.861136311594053, 0.347854845137454)
    }
};

LineGaussLegendreIntegrationPoints5::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.906179845938664, 0.236926885056189),
        IntegrationPointType(-0.538469310105683, 0.478628670499366),
        IntegrationPointType( 0.000000000000000, 0.568888888888889),
        IntegrationPointType( 0.538469310105683, 0.478628670499366),
        IntegrationPointType( 0.906179845938664, 0.236926885056189)
    }
};

LineGaussLegendreIntegrationPoints6::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints6::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.9324695142031521, 0.1713244923791704),
        IntegrationPointType(-0.6612093864662645, 0.3607615730481386),
        IntegrationPointType(-0.2386191860831969, 0.4679139345726910),
        IntegrationPointType( 0.2386191860831969, 0.4679139345726910),
        IntegrationPointType( 0.6612093864662645, 0.3607615730481386),
        IntegrationPointType( 0.9324695142031521, 0.1713244923791704)
    }
};

LineGaussLegendreIntegrationPoints7::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints7::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.9491079123427585, 0.1294849661688697),
        IntegrationPointType(-0.7415311855993945, 0.2797053914892766),
        IntegrationPointType(-0.4058451513773972, 0.3818300505051189),
        IntegrationPointType( 0.0000000000000000, 0.4179591836734694),
        IntegrationPointType( 0.4058451513773972, 0.3818300505051189),
        IntegrationPointType( 0.7415311855993945, 0.2797053914892766),
        IntegrationPointType( 0.9491079123427585, 0.1294849661688697)
    }
};

LineGaussLegendreIntegrationPoints8::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints8::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.9602898564975363, 0.1012285362903763),
        IntegrationPointType(-0.7966664774136267, 0.2223810344533745),
        IntegrationPointType(-0.5255324099163290, 0.3137066458778873),
        IntegrationPointType(-0.1834346424956498, 0.3626837833783620),
        IntegrationPointType( 0.1834346424956498, 0.3626837833783620),
        IntegrationPointType( 0.5255324099163290, 0.3137066458778873),
        IntegrationPointType( 0.7966664774136267, 0.2223810344533745),
        IntegrationPointType( 0.9602898564975363, 0.1012285362903763)
    }
};

LineGaussLegendreIntegrationPoints9::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints9::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.9681602395076261, 0.0812743883615744),
        IntegrationPointType(-0.8360311073266358, 0.1806481606948574),
        IntegrationPointType(-0.6133714327005904, 0.2606106964029354),
        IntegrationPointType(-0.3242534234038089, 0.3123470770400029),
        IntegrationPointType( 0.0000000000000000, 0.3302393550012598),
        IntegrationPointType( 0.3242534234038089, 0.3123470770400029),
        IntegrationPointType( 0.6133714327005904, 0.2606106964029354),
        IntegrationPointType( 0.8360311073266358, 0.1806481606948574),
        IntegrationPointType( 0.9681602395076261, 0.0812743883615744)
    }
};

LineGaussLegendreIntegrationPoints10::IntegrationPointsArrayType LineGaussLegendreIntegrationPoints10::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.9739065285171717, 0.0666713443086881),
        IntegrationPointType(-0.8650633666889845, 0.1494513491505806),
        IntegrationPointType(-0.6794095682990244, 0.2190863625159820),
        IntegrationPointType(-0.4333953941292472, 0.2692667193099963),
        IntegrationPointType(-0.1488743389816312, 0.2955242247147529),
        IntegrationPointType( 0.1488743389816312, 0.2955242247147529),
        IntegrationPointType( 0.4333953941292472, 0.2692667193099963),
        IntegrationPointType( 0.6794095682990244, 0.2190863625159820),
        IntegrationPointType( 0.8650633666889845, 0.1494513491505806),
        IntegrationPointType( 0.9739065285171717, 0.0666713443086881)
    }
};

//Gauss-Lobatto

LineGaussLobattoIntegrationPoints1::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType(0.00, 2.00)
    }
};

LineGaussLobattoIntegrationPoints2::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00, 1.00),
        IntegrationPointType( 1.00, 1.00)
    }
};

LineGaussLobattoIntegrationPoints3::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00, 1.00 / 3.00),
        IntegrationPointType( 0.00, 4.00 / 3.00),
        IntegrationPointType( 1.00, 1.00 / 3.00)
    }
};

LineGaussLobattoIntegrationPoints4::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00, 1.00 / 6.00),
        IntegrationPointType(-std::sqrt(5.00) / 5.00, 5.00 / 6.00),
        IntegrationPointType( std::sqrt(5.00) / 5.00, 5.00 / 6.00),
        IntegrationPointType( 1.00, 1.00 / 6.00)
    }
};

LineGaussLobattoIntegrationPoints5::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00, 0.10),
        IntegrationPointType(-std::sqrt(21.00) / 7.00, 49.00 / 90.00),
        IntegrationPointType( 0.00, 32.00 / 45.00),
        IntegrationPointType( std::sqrt(21.00) / 7.00, 49.00 / 90.00),
        IntegrationPointType( 1.00, 0.10)
    }
};

LineGaussLobattoIntegrationPoints6::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints6::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00, 1.00 / 15.00),
        IntegrationPointType(-std::sqrt((7.00+2.00*std::sqrt(7)) / 21.00), (14.00-std::sqrt(7)) / 30.00),
        IntegrationPointType(-std::sqrt((7.00-2.00*std::sqrt(7)) / 21.00), (14.00+std::sqrt(7)) / 30.0),
        IntegrationPointType( std::sqrt((7.00-2.00*std::sqrt(7)) / 21.00), (14.00+std::sqrt(7)) / 30.00),
        IntegrationPointType( std::sqrt((7.00+2.00*std::sqrt(7)) / 21.00), (14.00-std::sqrt(7)) / 30.00),
        IntegrationPointType( 1.00, 1.00 / 15.00)
    }
};

//TRIANGLE:

//Gauss-Legendre

TriangleGaussLegendreIntegrationPoints1::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 )
    }
};

TriangleGaussLegendreIntegrationPoints2::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints2::msIntegrationPoints =
{

    {
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 )
    }
};

TriangleGaussLegendreIntegrationPoints3::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType( 1.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
        IntegrationPointType( 3.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
        IntegrationPointType( 1.00 / 5.00 , 3.00 / 5.00 , 25.00 / 96.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -27.00 / 96.00 )
    }
};

#define KRATOS_TRI_G4_wa 0.109951743655322
#define KRATOS_TRI_G4_wb 0.223381589678011
#define KRATOS_TRI_G4_Na1 0.816847572980459
#define KRATOS_TRI_G4_Nb1 0.108103018168070
#define KRATOS_TRI_G4_Na2 0.091576213509771
#define KRATOS_TRI_G4_Nb2 0.445948490915965
TriangleGaussLegendreIntegrationPoints4::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType( KRATOS_TRI_G4_Na2, KRATOS_TRI_G4_Na2, KRATOS_TRI_G4_wa ),
        IntegrationPointType( KRATOS_TRI_G4_Na1, KRATOS_TRI_G4_Na2, KRATOS_TRI_G4_wa ),
        IntegrationPointType( KRATOS_TRI_G4_Na2, KRATOS_TRI_G4_Na1, KRATOS_TRI_G4_wa ),
        IntegrationPointType( KRATOS_TRI_G4_Nb2, KRATOS_TRI_G4_Nb2, KRATOS_TRI_G4_wb ),
        IntegrationPointType( KRATOS_TRI_G4_Nb1, KRATOS_TRI_G4_Nb2, KRATOS_TRI_G4_wb ),
        IntegrationPointType( KRATOS_TRI_G4_Nb2, KRATOS_TRI_G4_Nb1, KRATOS_TRI_G4_wb )
    }
};
#undef KRATOS_TRI_G4_wa
#undef KRATOS_TRI_G4_wb
#undef KRATOS_TRI_G4_Na1
#undef KRATOS_TRI_G4_Na2
#undef KRATOS_TRI_G4_Nb1
#undef KRATOS_TRI_G4_Nb2


//QUADRILATERAL:

//Gauss-Legendre

QuadrilateralGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.00 , 0.00 , 4.00 )
    }
};

QuadrilateralGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), 1.00 )
    }
};

QuadrilateralGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), 25.00/81.00 ),
        IntegrationPointType(             0.00 , -std::sqrt(3.00/5.00), 40.00/81.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), 25.00/81.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,             0.00, 40.00/81.00 ),
        IntegrationPointType(             0.00 ,             0.00, 64.00/81.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,             0.00, 40.00/81.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), 25.00/81.00 ),
        IntegrationPointType(             0.00 ,  std::sqrt(3.00/5.00), 40.00/81.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), 25.00/81.00 )
    }
};


QuadrilateralGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType( -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType( -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType( -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType( -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),

        IntegrationPointType( -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType( -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType( -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType( -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),

        IntegrationPointType(  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType(  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType(  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType(  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 + std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),

        IntegrationPointType(  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType(  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), -std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType(  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 - 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 + std::sqrt(5.0/6.0)/6.0)),
        IntegrationPointType(  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ),  std::sqrt( (3.0 + 2.0 * std::sqrt(6.0/5.0)) / 7.0 ), (0.5 - std::sqrt(5.0/6.0)/6.0)*(0.5 - std::sqrt(5.0/6.0)/6.0)),
    }
};


QuadrilateralGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
QuadrilateralGaussLegendreIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType( -0.906179845938664, -0.906179845938664, 0.236926885056189 * 0.236926885056189),
        IntegrationPointType( -0.906179845938664, -0.538469310105683, 0.236926885056189 * 0.478628670499366),
        IntegrationPointType( -0.906179845938664, 0.000000000000000, 0.236926885056189 * 0.568888888888889),
        IntegrationPointType( -0.906179845938664, 0.538469310105683, 0.236926885056189 * 0.478628670499366),
        IntegrationPointType( -0.906179845938664, 0.906179845938664, 0.236926885056189 * 0.236926885056189),

        IntegrationPointType( -0.538469310105683, -0.906179845938664, 0.478628670499366 * 0.236926885056189),
        IntegrationPointType( -0.538469310105683, -0.538469310105683, 0.478628670499366 * 0.478628670499366),
        IntegrationPointType( -0.538469310105683, 0.000000000000000, 0.478628670499366 * 0.568888888888889),
        IntegrationPointType( -0.538469310105683, 0.538469310105683, 0.478628670499366 * 0.478628670499366),
        IntegrationPointType( -0.538469310105683, 0.906179845938664, 0.478628670499366 * 0.236926885056189),

        IntegrationPointType( 0.000000000000000, -0.906179845938664, 0.568888888888889 * 0.236926885056189),
        IntegrationPointType( 0.000000000000000, -0.538469310105683, 0.568888888888889 * 0.478628670499366),
        IntegrationPointType( 0.000000000000000, 0.000000000000000, 0.568888888888889 * 0.568888888888889),
        IntegrationPointType( 0.000000000000000, 0.538469310105683, 0.568888888888889 * 0.478628670499366),
        IntegrationPointType( 0.000000000000000, 0.906179845938664, 0.568888888888889 * 0.236926885056189),

        IntegrationPointType( 0.538469310105683, -0.906179845938664, 0.478628670499366 * 0.236926885056189),
        IntegrationPointType( 0.538469310105683, -0.538469310105683, 0.478628670499366 * 0.478628670499366),
        IntegrationPointType( 0.538469310105683, 0.000000000000000, 0.478628670499366 * 0.568888888888889),
        IntegrationPointType( 0.538469310105683, 0.538469310105683, 0.478628670499366 * 0.478628670499366),
        IntegrationPointType( 0.538469310105683, 0.906179845938664, 0.478628670499366 * 0.236926885056189),

        IntegrationPointType( 0.906179845938664, -0.906179845938664, 0.236926885056189 * 0.236926885056189),
        IntegrationPointType( 0.906179845938664, -0.538469310105683, 0.236926885056189 * 0.478628670499366),
        IntegrationPointType( 0.906179845938664, 0.000000000000000, 0.236926885056189 * 0.568888888888889),
        IntegrationPointType( 0.906179845938664, 0.538469310105683, 0.236926885056189 * 0.478628670499366),
        IntegrationPointType( 0.906179845938664, 0.906179845938664, 0.236926885056189 * 0.236926885056189),
    }
};

//Gauss-Lobatto

QuadrilateralGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralGaussLobattoIntegrationPoints1::msIntegrationPoints =
{
	{
		IntegrationPointType( -1.00 , 0.00, 1.00 ),
		IntegrationPointType(  1.00 , 0.00, 1.00 )
	}
};

QuadrilateralGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralGaussLobattoIntegrationPoints2::msIntegrationPoints =
{
	{
		IntegrationPointType( -1.00 , -1.00, 0.50 ),
		IntegrationPointType(  1.00 , -1.00, 0.50 ),
		IntegrationPointType(  1.00 ,  1.00, 0.50 ),
		IntegrationPointType( -1.00 ,  1.00, 0.50 )
	}
};

//TETRAHEDRON:

//Gauss-Legendre

TetrahedronGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.25,0.25,0.25 , 1.00 / 6.00 )
    }
};

TetrahedronGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.58541020,0.13819660,0.13819660 , 1.00 / 24.00 ),
        IntegrationPointType( 0.13819660,0.58541020,0.13819660 , 1.00 / 24.00 ),
        IntegrationPointType( 0.13819660,0.13819660,0.58541020 , 1.00 / 24.00 ),
        IntegrationPointType( 0.13819660,0.13819660,0.13819660 , 1.00 / 24.00 )
    }
};

TetrahedronGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.25,0.25,0.25 , -0.1333333333333333333333333333333 ),
        IntegrationPointType( 0.5,0.1666666666666667,0.1666666666666667	, 0.075 ),
        IntegrationPointType( 0.1666666666666667,0.5,0.1666666666666667 , 0.075 ),
        IntegrationPointType( 0.1666666666666667,0.1666666666666667,0.5, 0.075 ),
        IntegrationPointType( 0.1666666666666667,0.1666666666666667,0.1666666666666667, 0.075 )
    }
};

// 	TetrahedronGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
// 			TetrahedronGaussLegendreIntegrationPoints4::msIntegrationPoints=
// 	{
// 		{
// 			IntegrationPointType( 1.0/4.0,  1.0/4.0,  1.0/4.0, -4.0/30.0 ),
// 			IntegrationPointType( 1.0/6.0,  1.0/6.0,  1.0/6.0, 9.0/120.0 ),
// 			IntegrationPointType( 1.0/2.0,  1.0/6.0,  1.0/6.0, 9.0/120.0 ),
// 			IntegrationPointType( 1.0/6.0,  1.0/2.0,  1.0/6.0, 9.0/120.0 ),
// 			IntegrationPointType( 1.0/6.0,  1.0/6.0,  1.0/2.0, 9.0/120.0 )
// 		}
// 	};


TetrahedronGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints4::msIntegrationPoints=
{
    {
        IntegrationPointType( tet10_a,  tet10_a,  tet10_a, -1.0/60.0 ),
        IntegrationPointType( tet10_c,  tet10_a,  tet10_a, -1.0/60.0 ),
        IntegrationPointType( tet10_a,  tet10_c,  tet10_a, -1.0/60.0 ),
        IntegrationPointType( tet10_a,  tet10_a,  tet10_c, -1.0/60.0 ),
        IntegrationPointType( tet10_b,  tet10_a,  tet10_a, -1.0/60.0 ),
        IntegrationPointType( tet10_b,  tet10_b,  tet10_a, -1.0/60.0 ),
        IntegrationPointType( tet10_a,  tet10_b,  tet10_a, -1.0/60.0 ),
        IntegrationPointType( tet10_a,  tet10_a,  tet10_b, -1.0/60.0 ),
        IntegrationPointType( tet10_b,  tet10_a,  tet10_b, -1.0/60.0 ),
        IntegrationPointType( tet10_a,  tet10_b,  tet10_b, -1.0/60.0 ),
    }
};


TetrahedronGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
TetrahedronGaussLegendreIntegrationPoints5::msIntegrationPoints=
{
    {
        IntegrationPointType(1.0/4.0, 1.0/4.0, 1.0/4.0,-74.0/5625.0 ),
        IntegrationPointType(1.0/14.0, 1.0/14.0, 1.0/14.0,343.0/45000.0 ),
        IntegrationPointType(11.0/14.0, 1.0/14.0, 1.0/14.0,343.0/45000.0 ),
        IntegrationPointType(1.0/14.0, 11.0/14.0, 1.0/14.0,343.0/45000.0 ),
        IntegrationPointType(1.0/14.0, 1.0/14.0, 11.0/14.0,343.0/45000.0 ),
        IntegrationPointType((1.0+std::sqrt(5.0/14.0))/4.0,(1.0-std::sqrt(5.0/14.0))/
        4.0, (1.0-std::sqrt(5.0/14.0))/4.0,56.0/2250.0 ),
        IntegrationPointType((1.0+std::sqrt(5.0/14.0))/4.0,(1.0+std::sqrt(5.0/14.0))/
        4.0, (1.0-std::sqrt(5.0/14.0))/4.0,56.0/2250.0 ),
        IntegrationPointType((1.0-std::sqrt(5.0/14.0))/4.0,(1.0+std::sqrt(5.0/14.0))/
        4.0, (1.0-std::sqrt(5.0/14.0))/4.0,56.0/2250.0 ),
        IntegrationPointType((1.0-std::sqrt(5.0/14.0))/4.0,(1.0-std::sqrt(5.0/14.0))/
        4.0, (1.0+std::sqrt(5.0/14.0))/4.0,56.0/2250.0 ),
        IntegrationPointType((1.0+std::sqrt(5.0/14.0))/4.0,(1.0-std::sqrt(5.0/14.0))/
        4.0, (1.0+std::sqrt(5.0/14.0))/4.0,56.0/2250.0 ),
        IntegrationPointType((1.0-std::sqrt(5.0/14.0))/4.0,(1.0+std::sqrt(5.0/14.0))/
        4.0, (1.0+std::sqrt(5.0/14.0))/4.0,56.0/2250.0 )
    }
};

//PRISM:

//Gauss-Legendre

// TODO: Check if the values are correct

PrismGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints1::msIntegrationPoints=
{
    {
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.50, 1.00 / 6.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.50, 1.00 / 6.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.50, 1.00 / 6.00)
    }
};

PrismGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints2::msIntegrationPoints=
{
    {
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, ((-1.0/std::sqrt(3.0)+1.0)/2.0), 1.00 / 12.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, ((-1.0/std::sqrt(3.0)+1.0)/2.0), 1.00 / 12.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, ((-1.0/std::sqrt(3.0)+1.0)/2.0), 1.00 / 12.00),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, (( 1.0/std::sqrt(3.0)+1.0)/2.0), 1.00 / 12.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, (( 1.0/std::sqrt(3.0)+1.0)/2.0), 1.00 / 12.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, (( 1.0/std::sqrt(3.0)+1.0)/2.0), 1.00 / 12.00)
    }
};

PrismGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints3::msIntegrationPoints=
{
    {
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, (-std::sqrt(3.00 / 5.00)+1.0)/2.0, 5.00 / 108.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, (-std::sqrt(3.00 / 5.00)+1.0)/2.0, 5.00 / 108.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, (-std::sqrt(3.00 / 5.00)+1.0)/2.0, 5.00 / 108.00),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.500000000000000000000000000000,  4.00 /  54.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.500000000000000000000000000000,  4.00 /  54.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.500000000000000000000000000000,  4.00 /  54.00),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, ( std::sqrt(3.00 / 5.00)+1.0)/2.0, 5.00 / 108.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, ( std::sqrt(3.00 / 5.00)+1.0)/2.0, 5.00 / 108.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, ( std::sqrt(3.00 / 5.00)+1.0)/2.0, 5.00 / 108.00)
    }
};

PrismGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints4::msIntegrationPoints=
{
    {
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.06943184420297371239, 0.02898790376145448812),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.06943184420297371239, 0.02898790376145448812),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.06943184420297371239, 0.02898790376145448812),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.33000947820757186760, 0.05434542957187884522),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.33000947820757186760, 0.05434542957187884522),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.33000947820757186760, 0.05434542957187884522),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.66999052179242813240, 0.05434542957187884522),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.66999052179242813240, 0.05434542957187884522),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.66999052179242813240, 0.05434542957187884522),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.93056815579702628761, 0.02898790376145448812),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.93056815579702628761, 0.02898790376145448812),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.93056815579702628761, 0.02898790376145448812)
    }
};

PrismGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints5::msIntegrationPoints=
{
    {
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.04691007703066800360, 0.01974390708801575729),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.04691007703066800360, 0.01974390708801575729),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.04691007703066800360, 0.01974390708801575729),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.23076534494715845448, 0.03988572254161387234),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.23076534494715845448, 0.03988572254161387234),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.23076534494715845448, 0.03988572254161387234),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.50000000000000000000, 0.04740740740740740741),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.50000000000000000000, 0.04740740740740740741),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.50000000000000000000, 0.04740740740740740741),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.76923465505284154552, 0.03988572254161387234),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.76923465505284154552, 0.03988572254161387234),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.76923465505284154552, 0.03988572254161387234),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.95308992296933199640, 0.01974390708801575729),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.95308992296933199640, 0.01974390708801575729),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.95308992296933199640, 0.01974390708801575729)
    }
};

//Gauss-Lobatto

PrismGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
PrismGaussLobattoIntegrationPoints1::msIntegrationPoints =
{
	{
		IntegrationPointType( 0.0 , 0.0 , 0.5 , 1.00 / 6.00 ),
		IntegrationPointType( 1.0 , 0.0 , 0.5 , 1.00 / 6.00 ),
		IntegrationPointType( 0.0 , 1.0 , 0.5 , 1.00 / 6.00 )
	}
};

PrismGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
PrismGaussLobattoIntegrationPoints2::msIntegrationPoints =
{
	{
		IntegrationPointType( 0.0 , 0.0 , 0.0 , 1.00 / 12.00 ),
		IntegrationPointType( 1.0 , 0.0 , 0.0 , 1.00 / 12.00 ),
		IntegrationPointType( 0.0 , 1.0 , 0.0 , 1.00 / 12.00 ),
		IntegrationPointType( 0.0 , 0.0 , 1.0 , 1.00 / 12.00 ),
		IntegrationPointType( 1.0 , 0.0 , 1.0 , 1.00 / 12.00 ),
		IntegrationPointType( 0.0 , 1.0 , 1.0 , 1.00 / 12.00 )
	}
};

//HEXAHEDRON:

//Gauss-Legendre

HexahedronGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.00 , 0.00, 0.00 , 8.00 )
    }
};

HexahedronGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 )
    }
};

HexahedronGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 ),
        IntegrationPointType(              0.0 , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 200.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,              0.0, -std::sqrt(3.00/5.00), 200.00/729.00 ),
        IntegrationPointType(              0.0 ,              0.0, -std::sqrt(3.00/5.00), 320.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,              0.0, -std::sqrt(3.00/5.00), 200.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 ),
        IntegrationPointType(              0.0 ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 200.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),              0.0, 200.00/729.00 ),
        IntegrationPointType(              0.0 , -std::sqrt(3.00/5.00),              0.0, 320.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),              0.0, 200.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,              0.0,              0.0, 320.00/729.00 ),
        IntegrationPointType(              0.0 ,              0.0,              0.0, 512.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,              0.0,              0.0, 320.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),              0.0, 200.00/729.00 ),
        IntegrationPointType(              0.0 ,  std::sqrt(3.00/5.00),              0.0, 320.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),              0.0, 200.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 ),
        IntegrationPointType(              0.0 , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 200.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,              0.0,  std::sqrt(3.00/5.00), 200.00/729.00 ),
        IntegrationPointType(              0.0 ,              0.0,  std::sqrt(3.00/5.00), 320.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,              0.0,  std::sqrt(3.00/5.00), 200.00/729.00 ),

        IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 ),
        IntegrationPointType(              0.0 ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 200.00/729.00 ),
        IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 )
    }
};

HexahedronGaussLegendreIntegrationPoints4::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454),
        IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454),
        IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454),
        IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454),
        IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990),
        IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454),
        IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454),
        IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087),
        IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454),
        IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098),
        IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454)
    }
};

HexahedronGaussLegendreIntegrationPoints5::IntegrationPointsArrayType
HexahedronGaussLegendreIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092),
        IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0 , -0.90617984593866399280 , -0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092),
        IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0 , -0.53846931010568309104 , -0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , 0 , -0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(-0.53846931010568309104 , 0 , -0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0 , 0 , -0.90617984593866399280 , 0.07667773006934522489),
        IntegrationPointType(0.53846931010568309104 , 0 , -0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0.90617984593866399280 , 0 , -0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0 , 0.53846931010568309104 , -0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092),
        IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0 , 0.90617984593866399280 , -0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092),
        IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0 , -0.90617984593866399280 , -0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0 , -0.53846931010568309104 , -0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.90617984593866399280 , 0 , -0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(-0.53846931010568309104 , 0 , -0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0 , 0 , -0.53846931010568309104 , 0.15490078296220484370),
        IntegrationPointType(0.53846931010568309104 , 0 , -0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0.90617984593866399280 , 0 , -0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0 , 0.53846931010568309104 , -0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0 , 0.90617984593866399280 , -0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , 0 , 0.031934207352848290676),
        IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , 0 , 0.06451200000000000000),
        IntegrationPointType(0 , -0.90617984593866399280 , 0 , 0.07667773006934522489),
        IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , 0 , 0.06451200000000000000),
        IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , 0 , 0.031934207352848290676),
        IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , 0 , 0.06451200000000000000),
        IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , 0 , 0.13032414106964827997),
        IntegrationPointType(0 , -0.53846931010568309104 , 0 , 0.15490078296220484370),
        IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , 0 , 0.13032414106964827997),
        IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , 0 , 0.06451200000000000000),
        IntegrationPointType(-0.90617984593866399280 , 0 , 0 , 0.07667773006934522489),
        IntegrationPointType(-0.53846931010568309104 , 0 , 0 , 0.15490078296220484370),
        IntegrationPointType(0 , 0 , 0 , 0.18411210973936899863),
        IntegrationPointType(0.53846931010568309104 , 0 , 0 , 0.15490078296220484370),
        IntegrationPointType(0.90617984593866399280 , 0 , 0 , 0.07667773006934522489),
        IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , 0 , 0.06451200000000000000),
        IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , 0 , 0.13032414106964827997),
        IntegrationPointType(0 , 0.53846931010568309104 , 0 , 0.15490078296220484370),
        IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , 0 , 0.13032414106964827997),
        IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , 0 , 0.06451200000000000000),
        IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , 0 , 0.031934207352848290676),
        IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , 0 , 0.06451200000000000000),
        IntegrationPointType(0 , 0.90617984593866399280 , 0 , 0.07667773006934522489),
        IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , 0 , 0.06451200000000000000),
        IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , 0 , 0.031934207352848290676),
        IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0 , -0.90617984593866399280 , 0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0 , -0.53846931010568309104 , 0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.90617984593866399280 , 0 , 0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(-0.53846931010568309104 , 0 , 0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0 , 0 , 0.53846931010568309104 , 0.15490078296220484370),
        IntegrationPointType(0.53846931010568309104 , 0 , 0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0.90617984593866399280 , 0 , 0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0 , 0.53846931010568309104 , 0.53846931010568309104 , 0.13032414106964827997),
        IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967),
        IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0 , 0.90617984593866399280 , 0.53846931010568309104 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092),
        IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0 , -0.90617984593866399280 , 0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092),
        IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0 , -0.53846931010568309104 , 0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , 0 , 0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(-0.53846931010568309104 , 0 , 0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0 , 0 , 0.90617984593866399280 , 0.07667773006934522489),
        IntegrationPointType(0.53846931010568309104 , 0 , 0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0.90617984593866399280 , 0 , 0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0 , 0.53846931010568309104 , 0.90617984593866399280 , 0.06451200000000000000),
        IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748),
        IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092),
        IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0 , 0.90617984593866399280 , 0.90617984593866399280 , 0.031934207352848290676),
        IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524),
        IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092)
    }
};

//Gauss-Lobatto

HexahedronGaussLobattoIntegrationPoints1::IntegrationPointsArrayType
HexahedronGaussLobattoIntegrationPoints1::msIntegrationPoints =
{
	{
		IntegrationPointType( -1.00 , -1.00, 0.00, 1.00 ),
		IntegrationPointType(  1.00 , -1.00, 0.00, 1.00 ),
		IntegrationPointType(  1.00 ,  1.00, 0.00, 1.00 ),
		IntegrationPointType( -1.00 ,  1.00, 0.00, 1.00 )
	}
};

HexahedronGaussLobattoIntegrationPoints2::IntegrationPointsArrayType
HexahedronGaussLobattoIntegrationPoints2::msIntegrationPoints =
{
	{
		IntegrationPointType( -1.00 , -1.00, -1.00, 0.50 ),
		IntegrationPointType(  1.00 , -1.00, -1.00, 0.50 ),
		IntegrationPointType(  1.00 ,  1.00, -1.00, 0.50 ),
		IntegrationPointType( -1.00 ,  1.00, -1.00, 0.50 ),
		IntegrationPointType( -1.00 , -1.00,  1.00, 0.50 ),
		IntegrationPointType(  1.00 , -1.00,  1.00, 0.50 ),
		IntegrationPointType(  1.00 ,  1.00,  1.00, 0.50 ),
		IntegrationPointType( -1.00 ,  1.00,  1.00, 0.50 )
	}
};


}
