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

LineGaussLobattoIntegrationPoints7::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints7::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00, 1.00 / 21.00),
        IntegrationPointType(-std::sqrt((5.00/11.00) + (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 - 7.00*std::sqrt(15.00)) / 350.00),
        IntegrationPointType(-std::sqrt((5.00/11.00) - (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 + 7.00*std::sqrt(15.00)) / 350.00),
        IntegrationPointType(0.00, 256.00/525.00),
        IntegrationPointType(std::sqrt((5.00/11.00) - (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 + 7.00*std::sqrt(15.00)) / 350.00),
        IntegrationPointType(std::sqrt((5.00/11.00) + (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 - 7.00*std::sqrt(15.00)) / 350.00),
        IntegrationPointType(-1.00, 1.00 / 21.00)
    }
};

LineGaussLobattoIntegrationPoints8::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints8::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00             , 0.035714285714286),
        IntegrationPointType(-0.871740148509607, 0.210704227143506),
        IntegrationPointType(-0.591700181433142, 0.341122692483504),
        IntegrationPointType(-0.209299217902479, 0.412458794658704),
        IntegrationPointType( 0.209299217902479, 0.412458794658704),
        IntegrationPointType( 0.591700181433142, 0.341122692483504),
        IntegrationPointType( 0.871740148509607, 0.210704227143506),
        IntegrationPointType( 1.00             , 0.035714285714286)
    }
};

LineGaussLobattoIntegrationPoints9::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints9::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00             , 0.027777777777778),
        IntegrationPointType(-0.899757995411460, 0.165495361560806),
        IntegrationPointType(-0.677186279510738, 0.274538712500162),
        IntegrationPointType(-0.363117463826178, 0.346428510973046),
        IntegrationPointType( 0.00             , 0.371519274376417),
        IntegrationPointType( 0.363117463826178, 0.346428510973046),
        IntegrationPointType( 0.677186279510738, 0.274538712500162),
        IntegrationPointType( 0.899757995411460, 0.165495361560806),
        IntegrationPointType( 1.00             , 0.027777777777778)
    }
};

LineGaussLobattoIntegrationPoints10::IntegrationPointsArrayType LineGaussLobattoIntegrationPoints10::msIntegrationPoints =
{
    {
        IntegrationPointType(-1.00             , 0.022222222222222),
        IntegrationPointType(-0.919533908166459, 0.133305990851070),
        IntegrationPointType(-0.738773865105505, 0.224889342063126),
        IntegrationPointType(-0.477924949810444, 0.292042683679684),
        IntegrationPointType(-0.165278957666387, 0.327539761183897),
        IntegrationPointType( 0.165278957666387, 0.327539761183897),
        IntegrationPointType( 0.477924949810444, 0.292042683679684),
        IntegrationPointType( 0.738773865105505, 0.224889342063126),
        IntegrationPointType( 0.919533908166459, 0.133305990851070),
        IntegrationPointType( 1.00             , 0.022222222222222)
    }
};

// Collocation
LineCollocationIntegrationPoints1::IntegrationPointsArrayType LineCollocationIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.666666666667,0.666666666667),
        IntegrationPointType(0,0.666666666667),
        IntegrationPointType(0.666666666667,0.666666666667)
    }
};

LineCollocationIntegrationPoints2::IntegrationPointsArrayType LineCollocationIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.800000000000,0.400000000000),
        IntegrationPointType(-0.400000000000,0.400000000000),
        IntegrationPointType(0.0000000000000,0.400000000000),
        IntegrationPointType(0.400000000000,0.400000000000),
        IntegrationPointType(0.800000000000,0.400000000000)
    }
};

LineCollocationIntegrationPoints3::IntegrationPointsArrayType LineCollocationIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.857142857143,0.285714285714),
        IntegrationPointType(-0.571428571429,0.285714285714),
        IntegrationPointType(-0.285714285714,0.285714285714),
        IntegrationPointType(0,0.285714285714),
        IntegrationPointType(0.285714285714,0.285714285714),
        IntegrationPointType(0.571428571429,0.285714285714),
        IntegrationPointType(0.857142857143,0.285714285714)
    }
};

LineCollocationIntegrationPoints4::IntegrationPointsArrayType LineCollocationIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.888888888889,0.222222222222),
        IntegrationPointType(-0.666666666667,0.222222222222),
        IntegrationPointType(-0.444444444444,0.222222222222),
        IntegrationPointType(-0.222222222222,0.222222222222),
        IntegrationPointType(0,0.222222222222),
        IntegrationPointType(0.222222222222,0.222222222222),
        IntegrationPointType(0.444444444444,0.222222222222),
        IntegrationPointType(0.666666666667,0.222222222222),
        IntegrationPointType(0.888888888889,0.222222222222)
    }
};

LineCollocationIntegrationPoints5::IntegrationPointsArrayType LineCollocationIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.909090909091,0.181818181818),
        IntegrationPointType(-0.727272727273,0.181818181818),
        IntegrationPointType(-0.545454545455,0.181818181818),
        IntegrationPointType(-0.363636363636,0.181818181818),
        IntegrationPointType(-0.181818181818,0.181818181818),
        IntegrationPointType(0,0.181818181818),
        IntegrationPointType(0.181818181818,0.181818181818),
        IntegrationPointType(0.363636363636,0.181818181818),
        IntegrationPointType(0.545454545455,0.181818181818),
        IntegrationPointType(0.727272727273,0.181818181818),
        IntegrationPointType(0.909090909091,0.181818181818)
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

#define KRATOS_TRI_G5_wa 0.025422453185103408460
#define KRATOS_TRI_G5_wb 0.058393137863189683013
#define KRATOS_TRI_G5_wc 0.041425537809186787597
#define KRATOS_TRI_G5_N1 0.87382197101699554332
#define KRATOS_TRI_G5_N2 0.063089014491502228340
#define KRATOS_TRI_G5_N3 0.50142650965817915742
#define KRATOS_TRI_G5_N4 0.24928674517091042129
#define KRATOS_TRI_G5_N5 0.053145049844816947353
#define KRATOS_TRI_G5_N6 0.31035245103378440542
#define KRATOS_TRI_G5_N7 0.63650249912139864723
TriangleGaussLegendreIntegrationPoints5::IntegrationPointsArrayType TriangleGaussLegendreIntegrationPoints5::msIntegrationPoints =
{
    {
	IntegrationPointType( KRATOS_TRI_G5_N1, KRATOS_TRI_G5_N2, KRATOS_TRI_G5_wa ),
        IntegrationPointType( KRATOS_TRI_G5_N2, KRATOS_TRI_G5_N1, KRATOS_TRI_G5_wa ),
        IntegrationPointType( KRATOS_TRI_G5_N2, KRATOS_TRI_G5_N2, KRATOS_TRI_G5_wa ),
        IntegrationPointType( KRATOS_TRI_G5_N3, KRATOS_TRI_G5_N4, KRATOS_TRI_G5_wb ),
        IntegrationPointType( KRATOS_TRI_G5_N4, KRATOS_TRI_G5_N3, KRATOS_TRI_G5_wb ),
        IntegrationPointType( KRATOS_TRI_G5_N4, KRATOS_TRI_G5_N4, KRATOS_TRI_G5_wb ),
	IntegrationPointType( KRATOS_TRI_G5_N5, KRATOS_TRI_G5_N6, KRATOS_TRI_G5_wc ),
        IntegrationPointType( KRATOS_TRI_G5_N6, KRATOS_TRI_G5_N5, KRATOS_TRI_G5_wc ),
        IntegrationPointType( KRATOS_TRI_G5_N5, KRATOS_TRI_G5_N7, KRATOS_TRI_G5_wc ),
        IntegrationPointType( KRATOS_TRI_G5_N6, KRATOS_TRI_G5_N7, KRATOS_TRI_G5_wc ),
        IntegrationPointType( KRATOS_TRI_G5_N7, KRATOS_TRI_G5_N5, KRATOS_TRI_G5_wc ),
        IntegrationPointType( KRATOS_TRI_G5_N7, KRATOS_TRI_G5_N6, KRATOS_TRI_G5_wc )
    }
};
#undef KRATOS_TRI_G5_wa
#undef KRATOS_TRI_G5_wb
#undef KRATOS_TRI_G5_wc
#undef KRATOS_TRI_G5_N1
#undef KRATOS_TRI_G5_N2
#undef KRATOS_TRI_G5_N3
#undef KRATOS_TRI_G5_N4
#undef KRATOS_TRI_G5_N5
#undef KRATOS_TRI_G5_N6
#undef KRATOS_TRI_G5_N7

TriangleGaussRadauIntegrationPoints1::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.5 , 0.5 , 1.00 / 3.00 ),
        IntegrationPointType( 0.5 , 0.0 , 1.00 / 3.00 ),
        IntegrationPointType( 0.0 , 0.5 , 1.00 / 3.00 )
    }
};

TriangleGaussRadauIntegrationPoints2::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -0.562500000000000 ),
        IntegrationPointType( 0.6 , 0.2 , 0.520833333333333 ),
        IntegrationPointType( 0.2 , 0.6 , 0.520833333333333 ),
        IntegrationPointType( 0.2 , 0.2 , 0.520833333333333 )
    }
};

TriangleGaussRadauIntegrationPoints3::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.816847572980459 , 0.091576213509771 , 0.109951743655322 ),
        IntegrationPointType( 0.091576213509771 , 0.816847572980459 , 0.109951743655322 ),
        IntegrationPointType( 0.091576213509771 , 0.091576213509771 , 0.109951743655322 ),
        IntegrationPointType( 0.108103018168070 , 0.445948490915965 , 0.223381589678011 ),
        IntegrationPointType( 0.445948490915965 , 0.108103018168070 , 0.223381589678011 ),
        IntegrationPointType( 0.445948490915965 , 0.445948490915965 , 0.223381589678011 )
    }
};

TriangleGaussRadauIntegrationPoints4::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 0.225000000000000 ),
        IntegrationPointType( 0.797426985353087 , 0.101286507323456 , 0.125939180544827 ),
        IntegrationPointType( 0.101286507323456 , 0.797426985353087 , 0.125939180544827 ),
        IntegrationPointType( 0.101286507323456 , 0.101286507323456 , 0.125939180544827 ),
        IntegrationPointType( 0.470142064105115 , 0.059715871789770 , 0.132394152788506 ),
        IntegrationPointType( 0.059715871789770 , 0.470142064105115 , 0.132394152788506 ),
        IntegrationPointType( 0.470142064105115 , 0.470142064105115 , 0.132394152788506 )
    }
};

TriangleGaussRadauIntegrationPoints5::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.873821971016996 , 0.063089014491502 , 0.050844906370207 ),
        IntegrationPointType( 0.063089014491502 , 0.873821971016996 , 0.050844906370207 ),
        IntegrationPointType( 0.063089014491502 , 0.063089014491502 , 0.050844906370207 ),
        IntegrationPointType( 0.501426509658179 , 0.249286745170910 , 0.116786275726379 ),
        IntegrationPointType( 0.249286745170910 , 0.501426509658179 , 0.116786275726379 ),
        IntegrationPointType( 0.249286745170910 , 0.249286745170910 , 0.116786275726379 ),
        IntegrationPointType( 0.636502499121399 , 0.310352451033785 , 0.082851075618374 ),
        IntegrationPointType( 0.636502499121399 , 0.053145049844816 , 0.082851075618374 ),
        IntegrationPointType( 0.310352451033785 , 0.636502499121399 , 0.082851075618374 ),
        IntegrationPointType( 0.310352451033785 , 0.053145049844816 , 0.082851075618374 ),
        IntegrationPointType( 0.053145049844816 , 0.636502499121399 , 0.082851075618374 ),
        IntegrationPointType( 0.053145049844816 , 0.310352451033785 , 0.082851075618374 )
    }
};

TriangleGaussRadauIntegrationPoints6::IntegrationPointsArrayType TriangleGaussRadauIntegrationPoints6::msIntegrationPoints =
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -0.149570044467670 ),
        IntegrationPointType( 0.479308067841923 , 0.260345966079038 , 0.175615257433204 ),
        IntegrationPointType( 0.260345966079038 , 0.479308067841923 , 0.175615257433204 ),
        IntegrationPointType( 0.260345966079038 , 0.260345966079038 , 0.175615257433204 ),
        IntegrationPointType( 0.869739794195568 , 0.065130102902216 , 0.053347235608839 ),
        IntegrationPointType( 0.065130102902216 , 0.869739794195568 , 0.053347235608839 ),
        IntegrationPointType( 0.065130102902216 , 0.065130102902216 , 0.053347235608839 ),
        IntegrationPointType( 0.638444188569809 , 0.312865496004875 , 0.077113760890257 ),
        IntegrationPointType( 0.638444188569809 , 0.048690315425316 , 0.077113760890257 ),
        IntegrationPointType( 0.312865496004875 , 0.638444188569809 , 0.077113760890257 ),
        IntegrationPointType( 0.312865496004875 , 0.048690315425316 , 0.077113760890257 ),
        IntegrationPointType( 0.048690315425316 , 0.638444188569809 , 0.077113760890257 ),
        IntegrationPointType( 0.048690315425316 , 0.312865496004875 , 0.077113760890257 )
    }
};

// Collocation
TriangleCollocationIntegrationPoints1::IntegrationPointsArrayType TriangleCollocationIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType(0.166666666667,0.166666666667,0.166666666667),
        IntegrationPointType(0.166666666667,0.666666666667,0.166666666667),
        IntegrationPointType(0.666666666667,0.666666666667,0.166666666667)
    }
};

TriangleCollocationIntegrationPoints2::IntegrationPointsArrayType TriangleCollocationIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType(0.111111111111,0.111111111111,0.0833333333333),
        IntegrationPointType(0.111111111111,0.444444444444,0.0833333333333),
        IntegrationPointType(0.111111111111,0.777777777778,0.0833333333333),
        IntegrationPointType(0.444444444444,0.444444444444,0.0833333333333),
        IntegrationPointType(0.444444444444,0.777777777778,0.0833333333333),
        IntegrationPointType(0.777777777778,0.777777777778,0.0833333333333)
    }
};

TriangleCollocationIntegrationPoints3::IntegrationPointsArrayType TriangleCollocationIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType(0.0833333333333,0.0833333333333,0.0500000000000),
        IntegrationPointType(0.0833333333333,0.333333333333,0.0500000000000),
        IntegrationPointType(0.0833333333333,0.583333333333,0.0500000000000),
        IntegrationPointType(0.0833333333333,0.833333333333,0.0500000000000),
        IntegrationPointType(0.333333333333,0.333333333333,0.0500000000000),
        IntegrationPointType(0.333333333333,0.583333333333,0.0500000000000),
        IntegrationPointType(0.333333333333,0.833333333333,0.0500000000000),
        IntegrationPointType(0.583333333333,0.583333333333,0.0500000000000),
        IntegrationPointType(0.583333333333,0.833333333333,0.0500000000000),
        IntegrationPointType(0.833333333333,0.833333333333,0.0500000000000)
    }
};

TriangleCollocationIntegrationPoints4::IntegrationPointsArrayType TriangleCollocationIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType(0.0666666666667,0.0666666666667,0.0333333333333),
        IntegrationPointType(0.0666666666667,0.266666666667,0.0333333333333),
        IntegrationPointType(0.0666666666667,0.466666666667,0.0333333333333),
        IntegrationPointType(0.0666666666667,0.666666666667,0.0333333333333),
        IntegrationPointType(0.0666666666667,0.866666666667,0.0333333333333),
        IntegrationPointType(0.266666666667,0.266666666667,0.0333333333333),
        IntegrationPointType(0.266666666667,0.466666666667,0.0333333333333),
        IntegrationPointType(0.266666666667,0.666666666667,0.0333333333333),
        IntegrationPointType(0.266666666667,0.866666666667,0.0333333333333),
        IntegrationPointType(0.466666666667,0.466666666667,0.0333333333333),
        IntegrationPointType(0.466666666667,0.666666666667,0.0333333333333),
        IntegrationPointType(0.466666666667,0.866666666667,0.0333333333333),
        IntegrationPointType(0.666666666667,0.666666666667,0.0333333333333),
        IntegrationPointType(0.666666666667,0.866666666667,0.0333333333333),
        IntegrationPointType(0.866666666667,0.866666666667,0.0333333333333)
    }
};

TriangleCollocationIntegrationPoints5::IntegrationPointsArrayType TriangleCollocationIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType(0.0555555555556,0.0555555555556,0.0238095238095),
        IntegrationPointType(0.0555555555556,0.222222222222,0.0238095238095),
        IntegrationPointType(0.0555555555556,0.388888888889,0.0238095238095),
        IntegrationPointType(0.0555555555556,0.555555555556,0.0238095238095),
        IntegrationPointType(0.0555555555556,0.722222222222,0.0238095238095),
        IntegrationPointType(0.0555555555556,0.888888888889,0.0238095238095),
        IntegrationPointType(0.222222222222,0.222222222222,0.0238095238095),
        IntegrationPointType(0.222222222222,0.388888888889,0.0238095238095),
        IntegrationPointType(0.222222222222,0.555555555556,0.0238095238095),
        IntegrationPointType(0.222222222222,0.722222222222,0.0238095238095),
        IntegrationPointType(0.222222222222,0.888888888889,0.0238095238095),
        IntegrationPointType(0.388888888889,0.388888888889,0.0238095238095),
        IntegrationPointType(0.388888888889,0.555555555556,0.0238095238095),
        IntegrationPointType(0.388888888889,0.722222222222,0.0238095238095),
        IntegrationPointType(0.388888888889,0.888888888889,0.0238095238095),
        IntegrationPointType(0.555555555556,0.555555555556,0.0238095238095),
        IntegrationPointType(0.555555555556,0.722222222222,0.0238095238095),
        IntegrationPointType(0.555555555556,0.888888888889,0.0238095238095),
        IntegrationPointType(0.722222222222,0.722222222222,0.0238095238095),
        IntegrationPointType(0.722222222222,0.888888888889,0.0238095238095),
        IntegrationPointType(0.888888888889,0.888888888889,0.0238095238095)
    }
};

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

// Collocation
QuadrilateralCollocationIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints1::msIntegrationPoints =
{
    {
IntegrationPointType(-0.500000000000,-0.500000000000,1.00000000000),
IntegrationPointType(-0.500000000000,0.500000000000,1.00000000000),
IntegrationPointType(0.500000000000,-0.500000000000,1.00000000000),
IntegrationPointType(0.500000000000,0.500000000000,1.00000000000)
    }
};

QuadrilateralCollocationIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints2::msIntegrationPoints =
{
    {
IntegrationPointType(-0.666666666667,-0.666666666667,0.444444444444),
IntegrationPointType(-0.666666666667,0,0.444444444444),
IntegrationPointType(-0.666666666667,0.666666666667,0.444444444444),
IntegrationPointType(0,-0.666666666667,0.444444444444),
IntegrationPointType(0,0,0.444444444444),
IntegrationPointType(0,0.666666666667,0.444444444444),
IntegrationPointType(0.666666666667,-0.666666666667,0.444444444444),
IntegrationPointType(0.666666666667,0,0.444444444444),
IntegrationPointType(0.666666666667,0.666666666667,0.444444444444)
    }
};

QuadrilateralCollocationIntegrationPoints3::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints3::msIntegrationPoints =
{
    {
IntegrationPointType(-0.750000000000,-0.750000000000,0.250000000000),
IntegrationPointType(-0.750000000000,-0.250000000000,0.250000000000),
IntegrationPointType(-0.750000000000,0.250000000000,0.250000000000),
IntegrationPointType(-0.750000000000,0.750000000000,0.250000000000),
IntegrationPointType(-0.250000000000,-0.750000000000,0.250000000000),
IntegrationPointType(-0.250000000000,-0.250000000000,0.250000000000),
IntegrationPointType(-0.250000000000,0.250000000000,0.250000000000),
IntegrationPointType(-0.250000000000,0.750000000000,0.250000000000),
IntegrationPointType(0.250000000000,-0.750000000000,0.250000000000),
IntegrationPointType(0.250000000000,-0.250000000000,0.250000000000),
IntegrationPointType(0.250000000000,0.250000000000,0.250000000000),
IntegrationPointType(0.250000000000,0.750000000000,0.250000000000),
IntegrationPointType(0.750000000000,-0.750000000000,0.250000000000),
IntegrationPointType(0.750000000000,-0.250000000000,0.250000000000),
IntegrationPointType(0.750000000000,0.250000000000,0.250000000000),
IntegrationPointType(0.750000000000,0.750000000000,0.250000000000)
    }
};


QuadrilateralCollocationIntegrationPoints4::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.800000000000,-0.800000000000,0.160000000000),
        IntegrationPointType(-0.800000000000,-0.400000000000,0.160000000000),
        IntegrationPointType(-0.800000000000, 0.000000000000,0.160000000000),
        IntegrationPointType(-0.800000000000, 0.400000000000,0.160000000000),
        IntegrationPointType(-0.800000000000, 0.800000000000,0.160000000000),
        IntegrationPointType(-0.400000000000,-0.800000000000,0.160000000000),
        IntegrationPointType(-0.400000000000,-0.400000000000,0.160000000000),
        IntegrationPointType(-0.400000000000, 0.000000000000,0.160000000000),
        IntegrationPointType(-0.400000000000, 0.400000000000,0.160000000000),
        IntegrationPointType(-0.400000000000, 0.800000000000,0.160000000000),
        IntegrationPointType( 0.000000000000,-0.800000000000,0.160000000000),
        IntegrationPointType( 0.000000000000,-0.400000000000,0.160000000000),
        IntegrationPointType( 0.000000000000, 0.000000000000,0.160000000000),
        IntegrationPointType( 0.000000000000, 0.400000000000,0.160000000000),
        IntegrationPointType( 0.000000000000, 0.800000000000,0.160000000000),
        IntegrationPointType(0.400000000000,-0.800000000000,0.160000000000),
        IntegrationPointType(0.400000000000,-0.400000000000,0.160000000000),
        IntegrationPointType(0.400000000000, 0.000000000000,0.160000000000),
        IntegrationPointType(0.400000000000,0.400000000000,0.160000000000),
        IntegrationPointType(0.400000000000,0.800000000000,0.160000000000),
        IntegrationPointType(0.800000000000,-0.800000000000,0.160000000000),
        IntegrationPointType(0.800000000000,-0.400000000000,0.160000000000),
        IntegrationPointType(0.800000000000, 0.000000000000, 0.160000000000),
        IntegrationPointType(0.800000000000,0.400000000000,0.160000000000),
        IntegrationPointType(0.800000000000,0.800000000000,0.160000000000)
    }
};


QuadrilateralCollocationIntegrationPoints5::IntegrationPointsArrayType
QuadrilateralCollocationIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.833333333333,-0.833333333333,0.111111111111),
        IntegrationPointType(-0.833333333333,-0.500000000000,0.111111111111),
        IntegrationPointType(-0.833333333333,-0.166666666667,0.111111111111),
        IntegrationPointType(-0.833333333333,0.166666666667,0.111111111111),
        IntegrationPointType(-0.833333333333,0.500000000000,0.111111111111),
        IntegrationPointType(-0.833333333333,0.833333333333,0.111111111111),
        IntegrationPointType(-0.500000000000,-0.833333333333,0.111111111111),
        IntegrationPointType(-0.500000000000,-0.500000000000,0.111111111111),
        IntegrationPointType(-0.500000000000,-0.166666666667,0.111111111111),
        IntegrationPointType(-0.500000000000,0.166666666667,0.111111111111),
        IntegrationPointType(-0.500000000000,0.500000000000,0.111111111111),
        IntegrationPointType(-0.500000000000,0.833333333333,0.111111111111),
        IntegrationPointType(-0.166666666667,-0.833333333333,0.111111111111),
        IntegrationPointType(-0.166666666667,-0.500000000000,0.111111111111),
        IntegrationPointType(-0.166666666667,-0.166666666667,0.111111111111),
        IntegrationPointType(-0.166666666667,0.166666666667,0.111111111111),
        IntegrationPointType(-0.166666666667,0.500000000000,0.111111111111),
        IntegrationPointType(-0.166666666667,0.833333333333,0.111111111111),
        IntegrationPointType(0.166666666667,-0.833333333333,0.111111111111),
        IntegrationPointType(0.166666666667,-0.500000000000,0.111111111111),
        IntegrationPointType(0.166666666667,-0.166666666667,0.111111111111),
        IntegrationPointType(0.166666666667,0.166666666667,0.111111111111),
        IntegrationPointType(0.166666666667,0.500000000000,0.111111111111),
        IntegrationPointType(0.166666666667,0.833333333333,0.111111111111),
        IntegrationPointType(0.500000000000,-0.833333333333,0.111111111111),
        IntegrationPointType(0.500000000000,-0.500000000000,0.111111111111),
        IntegrationPointType(0.500000000000,-0.166666666667,0.111111111111),
        IntegrationPointType(0.500000000000,0.166666666667,0.111111111111),
        IntegrationPointType(0.500000000000,0.500000000000,0.111111111111),
        IntegrationPointType(0.500000000000,0.833333333333,0.111111111111),
        IntegrationPointType(0.833333333333,-0.833333333333,0.111111111111),
        IntegrationPointType(0.833333333333,-0.500000000000,0.111111111111),
        IntegrationPointType(0.833333333333,-0.166666666667,0.111111111111),
        IntegrationPointType(0.833333333333,0.166666666667,0.111111111111),
        IntegrationPointType(0.833333333333,0.500000000000,0.111111111111),
        IntegrationPointType(0.833333333333,0.833333333333,0.111111111111)
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

/* EXTENDED VALUES (just one point along the plane) */

PrismGaussLegendreIntegrationPointsExt1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt1::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.2113248654051871177454, 1.00 / 4.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.7886751345948128822546, 1.00 / 4.00 )
    }
};

PrismGaussLegendreIntegrationPointsExt2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt2::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.1127016653792583114821, 5.00 / 36.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000, 4.00 / 18.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.8872983346207416885180, 5.00 / 36.00 )
    }
};

PrismGaussLegendreIntegrationPointsExt3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt3::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.0469100770306680036012, 0.0592317212640472718 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.2307653449471584544819, 0.1196571676248416170 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.5000000000000000000000, 0.1422222222222222222 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.7692346550528415455182, 0.1196571676248416170 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.9530899229693319963988, 0.0592317212640472718 )
    }
};

PrismGaussLegendreIntegrationPointsExt4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt4::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.0254460438286207377369, 0.0261224489795918367347 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.1292344072003027800681, 0.069926347872319166975  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.2970774243113014165467, 0.09545751262627973624   ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000, 0.1044897959183673469388 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.7029225756886985834533, 0.09545751262627973624   ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.8707655927996972199320, 0.069926347872319166975  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.9745539561713792622631, 0.0261224489795918367347 )
    }
};

PrismGaussLegendreIntegrationPointsExt5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsExt5::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.0108856709269715035981, 0.0139171417790434166207 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.0564687001159523504624, 0.031395092366226156159  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.1349239972129753379533, 0.0465725527319335628565 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.2404519353965940920372, 0.058298441147997619980  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.3652284220238275138343, 0.065701136127561665545  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000, 0.0682312716944751576786 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.6347715779761724861657, 0.065701136127561665545  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.7595480646034059079628, 0.058298441147997619980  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.8650760027870246620467, 0.0465725527319335628565 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.9435312998840476495376, 0.031395092366226156159  ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.9891143290730284964020, 0.0139171417790434166207 )
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
