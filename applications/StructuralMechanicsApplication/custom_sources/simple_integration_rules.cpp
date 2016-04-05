// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#include "custom_integration/simple_prism_gauss_legendre_integration_points.hpp"

namespace Kratos
{
    
//PRISM:

//Gauss-Legendre

/* SIMPLE VALUES (just one point along the plane) */

PrismGaussLegendreIntegrationPointsSimple1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple1::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.2113248654051871177454, 1.00 / 4.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.7886751345948128822546, 1.00 / 4.00 )
    }
};

PrismGaussLegendreIntegrationPointsSimple2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple2::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.1127016653792583114821, 5.00 / 36.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000, 4.00 / 18.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.8872983346207416885180, 5.00 / 36.00 )
    }
};

PrismGaussLegendreIntegrationPointsSimple3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple3::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.0469100770306680036012, 0.0592317212640472718 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.2307653449471584544819, 0.1196571676248416170 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.5000000000000000000000, 0.1422222222222222222 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.7692346550528415455182, 0.1196571676248416170 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.9530899229693319963988, 0.0592317212640472718 )
    }
};

PrismGaussLegendreIntegrationPointsSimple4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple4::msIntegrationPoints=
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

PrismGaussLegendreIntegrationPointsSimple5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple5::msIntegrationPoints=
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

}
