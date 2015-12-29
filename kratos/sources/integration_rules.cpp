//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   JMCarbonell $
//   Date:                $Date:   December 2015 $
//   Revision:            $Revision:         1.7 $
//
//

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

PrismGaussLegendreIntegrationPoints1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints1::msIntegrationPoints=
{
    {
        IntegrationPointType(0.25,0.25,0.5,1.0)
    }
};

PrismGaussLegendreIntegrationPoints2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints2::msIntegrationPoints=
{
    {
        IntegrationPointType(1.0/6.0,1.0/6.0,((-1.0/std::sqrt(3.0)+1.0)/2.0),1.0/6.0),
        IntegrationPointType(2.0/3.0,1.0/6.0,((-1.0/std::sqrt(3.0)+1.0)/2.0),1.0/6.0),
        IntegrationPointType(1.0/6.0,2.0/3.0,((-1.0/std::sqrt(3.0)+1.0)/2.0),1.0/6.0),
        IntegrationPointType(1.0/6.0,1.0/6.0,(( 1.0/std::sqrt(3.0)+1.0)/2.0),1.0/6.0),
        IntegrationPointType(2.0/3.0,1.0/6.0,(( 1.0/std::sqrt(3.0)+1.0)/2.0),1.0/6.0),
        IntegrationPointType(1.0/6.0,2.0/3.0,(( 1.0/std::sqrt(3.0)+1.0)/2.0),1.0/6.0)
    }
};

PrismGaussLegendreIntegrationPoints3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPoints3::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, -std::sqrt(3.00 / 5.00), 5.00 / 54.00),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, -std::sqrt(3.00 / 5.00), 5.00 / 54.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, -std::sqrt(3.00 / 5.00), 5.00 / 54.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, 0.0, 4.00 / 27.00),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, 0.0, 4.00 / 27.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, 0.0, 4.00 / 27.00),
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00, std::sqrt(3.00 / 5.00), 5.00 / 54.00),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00, std::sqrt(3.00 / 5.00), 5.00 / 54.00),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00, std::sqrt(3.00 / 5.00), 5.00 / 54.00)
    }
};

//Gauss-Lobatto

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

//Gauss-Lobatto

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
