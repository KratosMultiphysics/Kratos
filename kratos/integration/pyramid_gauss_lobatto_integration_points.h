//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{

class PyramidGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PyramidGaussLobattoIntegrationPoints1);

    static constexpr unsigned int Dimension = 3;

    typedef IntegrationPoint<Dimension> IntegrationPointType;

    typedef std::array<IntegrationPointType, 5> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 5;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( -1.0 , -1.0 , -1.0 , 5.0 / 12.0 ),
            IntegrationPointType(  1.0 , -1.0 , -1.0 , 5.0 / 12.0 ),
            IntegrationPointType(  1.0 ,  1.0 , -1.0 , 5.0 / 12.0 ),
            IntegrationPointType( -1.0 ,  1.0 , -1.0 , 5.0 / 12.0 ),
            IntegrationPointType(  0.0 ,  0.0 ,  1.0 , 1.0 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pyramid Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }

}; // Class PyramidGaussLobattoIntegrationPoints1

}



