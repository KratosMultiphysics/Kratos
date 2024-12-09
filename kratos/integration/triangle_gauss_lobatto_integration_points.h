//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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

class TriangleGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        const double one_over_six = 1.0 / 6.0;
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.0, 0.0 , one_over_six ),
            IntegrationPointType( 1.0, 0.0 , one_over_six ),
            IntegrationPointType( 0.0, 1.0 , one_over_six )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Lobatto quadrature 1 (3 points, degree of precision = 1)";
        return buffer.str();
    }


}; // Class TriangleGaussLobattoIntegrationPoints1

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.
