//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{

class TetrahedronGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedronGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        const double one_over_twenty_four = 1.0 / 24.0;
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.0, 0.0 , 0.0, one_over_twenty_four ),
            IntegrationPointType( 1.0, 0.0 , 0.0, one_over_twenty_four ),
            IntegrationPointType( 0.0, 1.0 , 0.0, one_over_twenty_four ),
            IntegrationPointType( 0.0, 0.0 , 1.0, one_over_twenty_four )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Tetrahedron Gauss-Lobatto quadrature 1 (4 points, degree of precision = 1)";
        return buffer.str();
    }
}; // Class TetrahedronGaussLobattoIntegrationPoints1

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.
