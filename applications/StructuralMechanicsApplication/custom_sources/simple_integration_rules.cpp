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

#include "custom_integration/simple_prism_gauss_legendre_integration_points.hpp"

namespace Kratos
{
    
//PRISM:

//Gauss-Legendre

/* SIMPLE VALUES (just one point along the plane) */

// TODO: Check if the values are correct

PrismGaussLegendreIntegrationPointsSimple1::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple1::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.5, 1.00 / 2.00 )
    }
};

PrismGaussLegendreIntegrationPointsSimple2::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple2::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.2113248654051871177454, 1.00 / 4.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.7886751345948128822546 , 1.00 / 4.00 )
    }
};

PrismGaussLegendreIntegrationPointsSimple3::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple3::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.1127016653792583114821, 5.00 / 36.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.5000000000000000000000, 4.00 / 18.00 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.8872983346207416885180, 5.00 / 36.00 )
    }
};

PrismGaussLegendreIntegrationPointsSimple4::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple4::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.0694318442029737123881, 0.0869637112843634643 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.3300094782075718675987, 0.1630362887156365357 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.6699905217924281324013, 0.1630362887156365357 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00, 0.9305681557970262876120, 0.0869637112843634643 )
    }
};

PrismGaussLegendreIntegrationPointsSimple5::IntegrationPointsArrayType
PrismGaussLegendreIntegrationPointsSimple5::msIntegrationPoints=
{
    {
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.0469100770306680036012, 0.0592317212640472718 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.2307653449471584544819, 0.1196571676248416170 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.5000000000000000000000, 0.1422222222222222222 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.7692346550528415455182, 0.1196571676248416170 ),
        IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00,  0.9530899229693319963988, 0.0592317212640472718 )
    }
};

}
