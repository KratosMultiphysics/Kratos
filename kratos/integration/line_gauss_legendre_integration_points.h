//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
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
class LineGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints1);
    
    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 1>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 1;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(0.00, 2.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints1


class LineGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints2);

    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 2>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 2;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-std::sqrt(1.00 / 3.00), 1.00),
            IntegrationPointType( std::sqrt(1.00 / 3.00), 1.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints2


class LineGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints3);
    
    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 3>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 3;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-std::sqrt(3.00 / 5.00), 5.00 / 9.00),
            IntegrationPointType( 0.00                  , 8.00 / 9.00),
            IntegrationPointType( std::sqrt(3.00 / 5.00), 5.00 / 9.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints3



class LineGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints4);

    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 4>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 4;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.861136311594053, 0.347854845137454),
            IntegrationPointType(-0.339981043584856, 0.652145154862546),
            IntegrationPointType( 0.339981043584856, 0.652145154862546),
            IntegrationPointType( 0.861136311594053, 0.347854845137454)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints4



class LineGaussLegendreIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints5);

    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 5>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 5;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.906179845938664, 0.236926885056189),
            IntegrationPointType(-0.538469310105683, 0.478628670499366),
            IntegrationPointType( 0.000000000000000, 0.568888888888889),
            IntegrationPointType( 0.538469310105683, 0.478628670499366),
            IntegrationPointType( 0.906179845938664, 0.236926885056189)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 5 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints4



class LineGaussLegendreIntegrationPoints6
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints6);
    
    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 6>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 6;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.9324695142031521, 0.1713244923791704),
            IntegrationPointType(-0.6612093864662645, 0.3607615730481386),
            IntegrationPointType(-0.2386191860831969, 0.4679139345726910),
            IntegrationPointType( 0.2386191860831969, 0.4679139345726910),
            IntegrationPointType( 0.6612093864662645, 0.3607615730481386),
            IntegrationPointType( 0.9324695142031521, 0.1713244923791704)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 6 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints6



class LineGaussLegendreIntegrationPoints7
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints7);
    
    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 7>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 7;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.9491079123427585, 0.1294849661688697),
            IntegrationPointType(-0.7415311855993945, 0.2797053914892766),
            IntegrationPointType(-0.4058451513773972, 0.3818300505051189),
            IntegrationPointType( 0.0000000000000000, 0.4179591836734694),
            IntegrationPointType( 0.4058451513773972, 0.3818300505051189),
            IntegrationPointType( 0.7415311855993945, 0.2797053914892766),
            IntegrationPointType( 0.9491079123427585, 0.1294849661688697)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 7 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints7



class LineGaussLegendreIntegrationPoints8
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints8);
    
    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 8>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 8;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.9602898564975363, 0.1012285362903763),
            IntegrationPointType(-0.7966664774136267, 0.2223810344533745),
            IntegrationPointType(-0.5255324099163290, 0.3137066458778873),
            IntegrationPointType(-0.1834346424956498, 0.3626837833783620),
            IntegrationPointType( 0.1834346424956498, 0.3626837833783620),
            IntegrationPointType( 0.5255324099163290, 0.3137066458778873),
            IntegrationPointType( 0.7966664774136267, 0.2223810344533745),
            IntegrationPointType( 0.9602898564975363, 0.1012285362903763)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 8 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints8



class LineGaussLegendreIntegrationPoints9
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints9);
    
    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 9>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 9;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.9681602395076261, 0.0812743883615744),
            IntegrationPointType(-0.8360311073266358, 0.1806481606948574),
            IntegrationPointType(-0.6133714327005904, 0.2606106964029354),
            IntegrationPointType(-0.3242534234038089, 0.3123470770400029),
            IntegrationPointType( 0.0000000000000000, 0.3302393550012598),
            IntegrationPointType( 0.3242534234038089, 0.3123470770400029),
            IntegrationPointType( 0.6133714327005904, 0.2606106964029354),
            IntegrationPointType( 0.8360311073266358, 0.1806481606948574),
            IntegrationPointType( 0.9681602395076261, 0.0812743883615744)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 9 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints9



class LineGaussLegendreIntegrationPoints10
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints10);
    
    static const unsigned int Dimension = 1;

    using IntegrationPointType = IntegrationPoint<1>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 10>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 10;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.9739065285171717, 0.0666713443086881),
            IntegrationPointType(-0.8650633666889845, 0.1494513491505806),
            IntegrationPointType(-0.6794095682990244, 0.2190863625159820),
            IntegrationPointType(-0.4333953941292472, 0.2692667193099963),
            IntegrationPointType(-0.1488743389816312, 0.2955242247147529),
            IntegrationPointType( 0.1488743389816312, 0.2955242247147529),
            IntegrationPointType( 0.4333953941292472, 0.2692667193099963),
            IntegrationPointType( 0.6794095682990244, 0.2190863625159820),
            IntegrationPointType( 0.8650633666889845, 0.1494513491505806),
            IntegrationPointType( 0.9739065285171717, 0.0666713443086881)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 10 ";
        return buffer.str();
    }


}; // Class LineGaussLegendreIntegrationPoints10

///@}


}  // namespace Kratos.


