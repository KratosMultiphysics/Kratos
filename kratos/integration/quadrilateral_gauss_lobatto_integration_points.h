//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Alejandro Cornejo
//
//

#pragma once


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

//TO BE COMPLETED: Only the needed ones have been implemented

namespace Kratos
{

class QuadrilateralGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 2> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 2;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( -1.00 , 0.00, 1.00 ),
            IntegrationPointType(  1.00 , 0.00, 1.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 1 ";
        return buffer.str();
    }


}; // Class QuadrilateralGaussLobattoIntegrationPoints1

class QuadrilateralGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( -1.00 , -1.00, 1.0 ),
            IntegrationPointType(  1.00 , -1.00, 1.0 ),
            IntegrationPointType(  1.00 ,  1.00, 1.0 ),
            IntegrationPointType( -1.00 ,  1.00, 1.0 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 2 ";
        return buffer.str();
    }

}; // Class QuadrilateralGaussLobattoIntegrationPoints2


class QuadrilateralGaussLobattoIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 9> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 9;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        // Define weights for 3rd-order Lobatto quadrature
        const double one_over_thirty_six  = 1.0 / 36.0;
        const double four_over_thirty_six = 4.0 * one_over_thirty_six;  // = 1/9
        const double sixteen_over_thirty_six = 16.0 * one_over_thirty_six;  // = 4/9

        // Define the integration points and their weights
        static const IntegrationPointsArrayType s_integration_points{{
            // Corner points
            IntegrationPointType( -1.00 , -1.00, one_over_thirty_six ),
            IntegrationPointType(  1.00 , -1.00, one_over_thirty_six ),
            IntegrationPointType( -1.00 ,  1.00, one_over_thirty_six ),
            IntegrationPointType(  1.00 ,  1.00, one_over_thirty_six ),

            // Edge midpoints
            IntegrationPointType(  0.00 , -1.00, four_over_thirty_six ),
            IntegrationPointType( -1.00 ,  0.00, four_over_thirty_six ),
            IntegrationPointType(  1.00 ,  0.00, four_over_thirty_six ),
            IntegrationPointType(  0.00 ,  1.00, four_over_thirty_six ),

            // Center point
            IntegrationPointType(  0.00 ,  0.00, sixteen_over_thirty_six )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 3 ";
        return buffer.str();
    }

}; // Class QuadrilateralGaussLobattoIntegrationPoints3

}


