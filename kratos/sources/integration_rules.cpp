/*
==============================================================================
Kratos
A General Purpose SoftKRATOS_TRI_G4_ware for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  softKRATOS_TRI_G4_ware  and  associated  documentation files  (the
"SoftKRATOS_TRI_G4_ware"), to  deal in  the SoftKRATOS_TRI_G4_ware without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  SoftKRATOS_TRI_G4_ware,  and to
permit persons to whom the SoftKRATOS_TRI_G4_ware  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the SoftKRATOS_TRI_G4_ware.

THE  SOFTKRATOS_TRI_G4_waRE IS  PROVIDED  "AS  IS", WITHOUT  KRATOS_TRI_G4_waRRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  KRATOS_TRI_G4_waRRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTKRATOS_TRI_G4_waRE OR THE USE OR OTHER DEALINGS IN THE SOFTKRATOS_TRI_G4_waRE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-10-23 11:35:13 $
//   Revision:            $Revision: 1.6 $
//
//


#include "integration/gauss_legendre_integration_points.h"
#include "integration/triangle_gaussian_integration_points.h"
#include "integration/tetrahedra_gaussian_integration_points.h"
#include "integration/hexahedra_gaussian_integration_points.h"
#include "integration/quadrilateral_gaussian_integration_points.h"
#include "integration/prism_gaussian_integration_points.h"
#include "integration/interface_integration_points.h"

#define tet10_a 0.108103018168070
#define tet10_b 0.445948490915965
#define tet10_c 0.816847572980459

namespace Kratos
{
GaussLegendreIntegrationPoints1::IntegrationPointsArrayType GaussLegendreIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType(0.00, 2.00)
    }
};

GaussLegendreIntegrationPoints2::IntegrationPointsArrayType GaussLegendreIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType(-std::sqrt(1.00 / 3.00), 1.00),
        IntegrationPointType( std::sqrt(1.00 / 3.00), 1.00)
    }
};

GaussLegendreIntegrationPoints3::IntegrationPointsArrayType GaussLegendreIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType(-std::sqrt(3.00 / 5.00), 5.00 / 9.00),
        IntegrationPointType( 0.00                  , 8.00 / 9.00),
        IntegrationPointType( std::sqrt(3.00 / 5.00), 5.00 / 9.00)
    }
};

GaussLegendreIntegrationPoints4::IntegrationPointsArrayType GaussLegendreIntegrationPoints4::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.861136311594053, 0.347854845137454),
        IntegrationPointType(-0.339981043584856, 0.652145154862546),
        IntegrationPointType( 0.339981043584856, 0.652145154862546),
        IntegrationPointType( 0.861136311594053, 0.347854845137454)
    }
};

GaussLegendreIntegrationPoints5::IntegrationPointsArrayType GaussLegendreIntegrationPoints5::msIntegrationPoints =
{
    {
        IntegrationPointType(-0.906179845938664, 0.236926885056189),
        IntegrationPointType(-0.538469310105683, 0.478628670499366),
        IntegrationPointType( 0.000000000000000, 0.568888888888889),
        IntegrationPointType( 0.538469310105683, 0.478628670499366),
        IntegrationPointType( 0.906179845938664, 0.236926885056189)
    }
};

TriangleGaussianIntegrationPoints1::IntegrationPointsArrayType TriangleGaussianIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 )
    }
};

TriangleGaussianIntegrationPoints2::IntegrationPointsArrayType TriangleGaussianIntegrationPoints2::msIntegrationPoints =
{

    {
        IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
        IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
        IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 )
    }
};

TriangleGaussianIntegrationPoints3::IntegrationPointsArrayType TriangleGaussianIntegrationPoints3::msIntegrationPoints =
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
TriangleGaussianIntegrationPoints4::IntegrationPointsArrayType TriangleGaussianIntegrationPoints4::msIntegrationPoints =
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

//tetrahedra 1 GP
TetrahedraGaussianIntegrationPoints1::IntegrationPointsArrayType
TetrahedraGaussianIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.25,0.25,0.25 , 1.00 / 6.00 )
    }
};

//tetrahedra 4 GP
TetrahedraGaussianIntegrationPoints2::IntegrationPointsArrayType
TetrahedraGaussianIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.58541020,0.13819660,0.13819660 , 1.00 / 24.00 ),
        IntegrationPointType( 0.13819660,0.58541020,0.13819660 , 1.00 / 24.00 ),
        IntegrationPointType( 0.13819660,0.13819660,0.58541020 , 1.00 / 24.00 ),
        IntegrationPointType( 0.13819660,0.13819660,0.13819660 , 1.00 / 24.00 )
    }
};

//tetrahedra 5 GP
TetrahedraGaussianIntegrationPoints3::IntegrationPointsArrayType
TetrahedraGaussianIntegrationPoints3::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.25,0.25,0.25 , -0.1333333333333333333333333333333 ),
        IntegrationPointType( 0.5,0.1666666666666667,0.1666666666666667	, 0.075 ),
        IntegrationPointType( 0.1666666666666667,0.5,0.1666666666666667 , 0.075 ),
        IntegrationPointType( 0.1666666666666667,0.1666666666666667,0.5, 0.075 ),
        IntegrationPointType( 0.1666666666666667,0.1666666666666667,0.1666666666666667, 0.075 )
    }
};

// 	TetrahedraGaussianIntegrationPoints4::IntegrationPointsArrayType
// 			TetrahedraGaussianIntegrationPoints4::msIntegrationPoints=
// 	{
// 		{
// 			IntegrationPointType( 1.0/4.0,  1.0/4.0,  1.0/4.0, -4.0/30.0 ),
// 			IntegrationPointType( 1.0/6.0,  1.0/6.0,  1.0/6.0, 9.0/120.0 ),
// 			IntegrationPointType( 1.0/2.0,  1.0/6.0,  1.0/6.0, 9.0/120.0 ),
// 			IntegrationPointType( 1.0/6.0,  1.0/2.0,  1.0/6.0, 9.0/120.0 ),
// 			IntegrationPointType( 1.0/6.0,  1.0/6.0,  1.0/2.0, 9.0/120.0 )
// 		}
// 	};

//tetrahedra 10 GP
TetrahedraGaussianIntegrationPoints4::IntegrationPointsArrayType
TetrahedraGaussianIntegrationPoints4::msIntegrationPoints=
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

//tetrahedra 11 GP
TetrahedraGaussianIntegrationPoints5::IntegrationPointsArrayType
TetrahedraGaussianIntegrationPoints5::msIntegrationPoints=
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

//Prism
PrismGaussianIntegrationPoints1::IntegrationPointsArrayType
PrismGaussianIntegrationPoints1::msIntegrationPoints=
{
    {
        IntegrationPointType(0.25,0.25,0.5,1.0)
    }
};

PrismGaussianIntegrationPoints2::IntegrationPointsArrayType
PrismGaussianIntegrationPoints2::msIntegrationPoints=
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

PrismGaussianIntegrationPoints3::IntegrationPointsArrayType
PrismGaussianIntegrationPoints3::msIntegrationPoints=
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


//quadrilateral!
QuadrilateralGaussianIntegrationPoints1::IntegrationPointsArrayType
QuadrilateralGaussianIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.00 , 0.00 , 4.00 )
    }
};

QuadrilateralGaussianIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralGaussianIntegrationPoints2::msIntegrationPoints =
{
    {
        IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), 1.00 ),
        IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), 1.00 )
    }
};

QuadrilateralGaussianIntegrationPoints3::IntegrationPointsArrayType
QuadrilateralGaussianIntegrationPoints3::msIntegrationPoints =
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


QuadrilateralGaussianIntegrationPoints4::IntegrationPointsArrayType
QuadrilateralGaussianIntegrationPoints4::msIntegrationPoints =
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


//hexahedra!
HexahedraGaussianIntegrationPoints1::IntegrationPointsArrayType
HexahedraGaussianIntegrationPoints1::msIntegrationPoints =
{
    {
        IntegrationPointType( 0.00 , 0.00, 0.00 , 8.00 )
    }
};

HexahedraGaussianIntegrationPoints2::IntegrationPointsArrayType
HexahedraGaussianIntegrationPoints2::msIntegrationPoints =
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

HexahedraGaussianIntegrationPoints3::IntegrationPointsArrayType
HexahedraGaussianIntegrationPoints3::msIntegrationPoints =
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

// interface geometries
QuadrilateralInterfaceLobattoIntegrationPoints2::IntegrationPointsArrayType
QuadrilateralInterfaceLobattoIntegrationPoints2::msIntegrationPoints =
{
	{
		IntegrationPointType( -1.00 , -1.00, 0.50 ),
		IntegrationPointType(  1.00 , -1.00, 0.50 ),
		IntegrationPointType(  1.00 ,  1.00, 0.50 ),
		IntegrationPointType( -1.00 ,  1.00, 0.50 )
	}
};

HexaedralInterfaceLobattoIntegrationPoints2::IntegrationPointsArrayType
HexaedralInterfaceLobattoIntegrationPoints2::msIntegrationPoints =
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

PrismInterfaceLobattoIntegrationPoints2::IntegrationPointsArrayType
PrismInterfaceLobattoIntegrationPoints2::msIntegrationPoints =
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

}