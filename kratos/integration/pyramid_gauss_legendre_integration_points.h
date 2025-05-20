//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Ashish Darekar
//
//  Reference:       (https://github.com/nschloe/quadpy/blob/main/src/quadpy/p3/_felippa.py)

#pragma once

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos {

class PyramidGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PyramidGaussLegendreIntegrationPoints1);

    static const unsigned int Dimension = 3;

    using IntegrationPointType = IntegrationPoint<3>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 1>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 1;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.0 , 0.0 , -0.5 , 4.7407407407407405 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pyramid Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }

}; // Class PyramidGaussLegendreIntegrationPoints1

class PyramidGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PyramidGaussLegendreIntegrationPoints2);
    
    static const unsigned int Dimension = 3;

    using IntegrationPointType = IntegrationPoint<3>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 5>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 5;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.5842373946721772 , 0.5842373946721772 , -0.6666666666666666 , 0.81 ),
            IntegrationPointType( -0.5842373946721772 , 0.5842373946721772 , -0.6666666666666666 , 0.81 ),
            IntegrationPointType( 0.5842373946721772 , -0.5842373946721772 , -0.6666666666666666 , 0.81 ),
            IntegrationPointType( -0.5842373946721772 , -0.5842373946721772 , -0.6666666666666666 , 0.81 ),
            IntegrationPointType( 0.0 , 0.0 , 0.4 , 4.62962962962963 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pyramid Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }

}; // Class PyramidGaussLegendreIntegrationPoints2

class PyramidGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PyramidGaussLegendreIntegrationPoints3);
    
    static const unsigned int Dimension = 3;

    using IntegrationPointType = IntegrationPoint<3>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 8>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 8;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.5773502691896257 , 0.5773502691896257 , 0.08830368802245057 , 0.9700392164449294 ),
            IntegrationPointType( -0.5773502691896257 , 0.5773502691896257 , 0.08830368802245057 , 0.9700392164449294 ),
            IntegrationPointType( 0.5773502691896257 , -0.5773502691896257 , 0.08830368802245057 , 0.9700392164449294 ),
            IntegrationPointType( -0.5773502691896257 , -0.5773502691896257 , 0.08830368802245057 , 0.9700392164449294 ),
            IntegrationPointType( 0.5773502691896257 , 0.5773502691896257 , -0.7549703546891172 , 0.6040348576291448 ),
            IntegrationPointType( -0.5773502691896257 , 0.5773502691896257 , -0.7549703546891172 , 0.6040348576291448 ),
            IntegrationPointType( 0.5773502691896257 , -0.5773502691896257 , -0.7549703546891172 , 0.6040348576291448 ),
            IntegrationPointType( -0.5773502691896257 , -0.5773502691896257 , -0.7549703546891172 , 0.6040348576291448 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pyramid Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }

}; // Class PyramidGaussLegendreIntegrationPoints3

class PyramidGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PyramidGaussLegendreIntegrationPoints4);
    
    static const unsigned int Dimension = 3;

    using IntegrationPointType = IntegrationPoint<3>;

    using IntegrationPointsArrayType = std::array<IntegrationPointType, 18>;

    using PointType = IntegrationPointType::PointType;

    static std::size_t IntegrationPointsNumber()
    {
        return 18;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.7745966692414834 , 0.7745966692414834 , 0.08830368802245057 , 0.2993948198904103 ),
            IntegrationPointType( -0.7745966692414834 , 0.7745966692414834 , 0.08830368802245057 , 0.2993948198904103 ),
            IntegrationPointType( 0.7745966692414834 , -0.7745966692414834 , 0.08830368802245057 , 0.2993948198904103 ),
            IntegrationPointType( -0.7745966692414834 , -0.7745966692414834 , 0.08830368802245057 , 0.2993948198904103 ),
            IntegrationPointType( 0.7745966692414834 , 0.0 , 0.08830368802245057 , 0.4790317118246565 ),
            IntegrationPointType( -0.7745966692414834 , 0.0 , 0.08830368802245057 , 0.4790317118246565 ),
            IntegrationPointType( 0.0 , 0.7745966692414834 , 0.08830368802245057 , 0.4790317118246565 ),
            IntegrationPointType( 0.0 , -0.7745966692414834 , 0.08830368802245057 , 0.4790317118246565 ),
            IntegrationPointType( 0.0 , 0.0 , 0.08830368802245057 , 0.7664507389194504 ),
            IntegrationPointType( 0.7745966692414834 , 0.7745966692414834 , -0.7549703546891172 , 0.18643051161393356 ),
            IntegrationPointType( -0.7745966692414834 , 0.7745966692414834 , -0.7549703546891172 , 0.18643051161393356 ),
            IntegrationPointType( 0.7745966692414834 , -0.7745966692414834 , -0.7549703546891172 , 0.18643051161393356 ),
            IntegrationPointType( -0.7745966692414834 , -0.7745966692414834 , -0.7549703546891172 , 0.18643051161393356 ),
            IntegrationPointType( 0.7745966692414834 , 0.0 , -0.7549703546891172 , 0.29828881858229367 ),
            IntegrationPointType( -0.7745966692414834 , 0.0 , -0.7549703546891172 , 0.29828881858229367 ),
            IntegrationPointType( 0.0 , 0.7745966692414834 , -0.7549703546891172 , 0.29828881858229367 ),
            IntegrationPointType( 0.0 , -0.7745966692414834 , -0.7549703546891172 , 0.29828881858229367 ),
            IntegrationPointType( 0.0 , 0.0 , -0.7549703546891172 , 0.4772621097316699 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pyramid Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }

}; // Class PyramidGaussLegendreIntegrationPoints4

class PyramidGaussLegendreIntegrationPoints5
{
public:
   KRATOS_CLASS_POINTER_DEFINITION(PyramidGaussLegendreIntegrationPoints5);
   
   static const unsigned int Dimension = 3;

   using IntegrationPointType = IntegrationPoint<3>;

   using IntegrationPointsArrayType = std::array<IntegrationPointType, 27>;

   using PointType = IntegrationPointType::PointType;

   static std::size_t IntegrationPointsNumber()
   {
       return 27;
   }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.7745966692414834 , 0.7745966692414834 , -0.8540119518537005 , 0.11287470009663707 ),
            IntegrationPointType( -0.7745966692414834 , 0.7745966692414834 , -0.8540119518537005 , 0.11287470009663707 ),
            IntegrationPointType( 0.7745966692414834 , -0.7745966692414834 , -0.8540119518537005 , 0.11287470009663707 ),
            IntegrationPointType( -0.7745966692414834 , -0.7745966692414834 , -0.8540119518537005 , 0.11287470009663707 ),
            IntegrationPointType( 0.7745966692414834 , 0.0 , -0.8540119518537005 , 0.1805995201546193 ),
            IntegrationPointType( -0.7745966692414834 , 0.0 , -0.8540119518537005 , 0.1805995201546193 ),
            IntegrationPointType( 0.0 , 0.7745966692414834 , -0.8540119518537005 , 0.1805995201546193 ),
            IntegrationPointType( 0.0 , -0.7745966692414834 , -0.8540119518537005 , 0.1805995201546193 ),
            IntegrationPointType( 0.0 , 0.0 , -0.8540119518537005 , 0.2889592322473909 ),
            IntegrationPointType( 0.7745966692414834 , 0.7745966692414834 , -0.3059924679232962 , 0.211713439795844 ),
            IntegrationPointType( -0.7745966692414834 , 0.7745966692414834 , -0.3059924679232962 , 0.211713439795844 ),
            IntegrationPointType( 0.7745966692414834 , -0.7745966692414834 , -0.3059924679232962 , 0.211713439795844 ),
            IntegrationPointType( -0.7745966692414834 , -0.7745966692414834 , -0.3059924679232962 , 0.211713439795844 ),
            IntegrationPointType( 0.7745966692414834 , 0.0 , -0.3059924679232962 , 0.33874150367335043 ),
            IntegrationPointType( -0.7745966692414834 , 0.0 , -0.3059924679232962 , 0.33874150367335043 ),
            IntegrationPointType( 0.0 , 0.7745966692414834 , -0.3059924679232962 , 0.33874150367335043 ),
            IntegrationPointType( 0.0 , -0.7745966692414834 , -0.3059924679232962 , 0.33874150367335043 ),
            IntegrationPointType( 0.0 , 0.0 , -0.3059924679232962 , 0.5419864058773607 ),
            IntegrationPointType( 0.7745966692414834 , 0.7745966692414834 , 0.41000441977699675 , 0.21244889714455584 ),
            IntegrationPointType( -0.7745966692414834 , 0.7745966692414834 , 0.41000441977699675 , 0.21244889714455584 ),
            IntegrationPointType( 0.7745966692414834 , -0.7745966692414834 , 0.41000441977699675 , 0.21244889714455584 ),
            IntegrationPointType( -0.7745966692414834 , -0.7745966692414834 , 0.41000441977699675 , 0.21244889714455584 ),
            IntegrationPointType( 0.7745966692414834 , 0.0 , 0.41000441977699675 , 0.33991823543128935 ),
            IntegrationPointType( -0.7745966692414834 , 0.0 , 0.41000441977699675 , 0.33991823543128935 ),
            IntegrationPointType( 0.0 , 0.7745966692414834 , 0.41000441977699675 , 0.33991823543128935 ),
            IntegrationPointType( 0.0 , -0.7745966692414834 , 0.41000441977699675 , 0.33991823543128935 ),
            IntegrationPointType( 0.0 , 0.0 , 0.41000441977699675 , 0.543869176690063 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pyramid Gauss-Legendre quadrature 5 ";
        return buffer.str();
    }

}; // Class PyramidGaussLegendreIntegrationPoints5

}  // namespace Kratos.
