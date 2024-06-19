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
        const double one_third = 1.0 / 3.0;
        const double four_thirds = 4.0 * one_third;
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( -1.00 ,  1.00, one_third ),
            IntegrationPointType( -1.00 , -1.00, one_third ),
            IntegrationPointType(  1.00 , -1.00, one_third ),
            IntegrationPointType(  1.00 ,  1.00, one_third ),

            IntegrationPointType(  0.00 ,  1.00, four_thirds ),
            IntegrationPointType( -1.00 ,  0.00, four_thirds ),
            IntegrationPointType(  0.00 , -1.00, four_thirds ),
            IntegrationPointType(  1.00 ,  0.00, four_thirds ),

            IntegrationPointType(  0.00 , 0.00, four_thirds )
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


