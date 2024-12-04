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
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{

class TriangleGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 1> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }


}; // Class TriangleGaussLegendreIntegrationPoints1


class TriangleGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints2);
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
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
            IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
            IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }


}; // Class TriangleGaussLegendreIntegrationPoints2


class TriangleGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLobattoIntegrationPoints2);
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
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.0, 0.0 , 1.00 / 6.00 ),
            IntegrationPointType( 1.0, 0.0 , 1.00 / 6.00 ),
            IntegrationPointType( 0.0, 1.0 , 1.00 / 6.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Lobatto quadrature 2 ";
        return buffer.str();
    }


}; // Class TriangleGaussLobattoIntegrationPoints2

class TriangleGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints3);
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
            IntegrationPointType( 1.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
            IntegrationPointType( 3.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 ),
            IntegrationPointType( 1.00 / 5.00 , 3.00 / 5.00 , 25.00 / 96.00 ),
            IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -27.00 / 96.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }


}; // Class TriangleGaussLegendreIntegrationPoints2


class TriangleGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        const double wa = 0.054975871827661;
        const double wb = 0.1116907948390055;
        const double Na1 = 0.816847572980459;
        const double Nb1 = 0.108103018168070;
        const double Na2 = 0.091576213509771;
        const double Nb2 = 0.445948490915965;

        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( Na2, Na2, wa ),
            IntegrationPointType( Na1, Na2, wa ),
            IntegrationPointType( Na2, Na1, wa ),
            IntegrationPointType( Nb2, Nb2, wb ),
            IntegrationPointType( Nb1, Nb2, wb ),
            IntegrationPointType( Nb2, Nb1, wb )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }


}; // Class TriangleGaussLegendreIntegrationPoints4

class TriangleGaussLegendreIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 12> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 12;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        const double wa = 0.025422453185103408460;
        const double wb = 0.058393137863189683013;
        const double wc = 0.041425537809186787597;

        const double N1 = 0.87382197101699554332;
        const double N2 = 0.063089014491502228340;
        const double N3 = 0.50142650965817915742;
        const double N4 = 0.24928674517091042129;
        const double N5 = 0.053145049844816947353;
        const double N6 = 0.31035245103378440542;
        const double N7 = 0.63650249912139864723;

        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( N1, N2, wa ),
            IntegrationPointType( N2, N1, wa ),
            IntegrationPointType( N2, N2, wa ),
            IntegrationPointType( N3, N4, wb ),
            IntegrationPointType( N4, N3, wb ),
            IntegrationPointType( N4, N4, wb ),
            IntegrationPointType( N5, N6, wc ),
            IntegrationPointType( N6, N5, wc ),
            IntegrationPointType( N5, N7, wc ),
            IntegrationPointType( N6, N7, wc ),
            IntegrationPointType( N7, N5, wc ),
            IntegrationPointType( N7, N6, wc )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 5 ";
        return buffer.str();
    }


}; // Class TriangleGaussLegendreIntegrationPoints5


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

