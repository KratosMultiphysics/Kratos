//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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

class LineGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints1);

    static constexpr unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<Dimension>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 2>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 2;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.0, 1.0),
            IntegrationPointType( 1.0, 1.0)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }

}; // Class LineGaussLobattoIntegrationPoints1

}



