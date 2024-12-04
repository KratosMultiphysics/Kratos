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

#pragma once


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{

class QuadrilateralGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLegendreIntegrationPoints1);
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
            IntegrationPointType( 0.00 , 0.00 , 4.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }


}; // Class QuadrilateralGaussLegendreIntegrationPoints1

class QuadrilateralGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLegendreIntegrationPoints2);
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
        const double one_over_sqrt_3 = 1.00 / std::sqrt(3.0);
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( -one_over_sqrt_3 , -one_over_sqrt_3, 1.00 ),
            IntegrationPointType(  one_over_sqrt_3 , -one_over_sqrt_3, 1.00 ),
            IntegrationPointType(  one_over_sqrt_3 ,  one_over_sqrt_3, 1.00 ),
            IntegrationPointType( -one_over_sqrt_3 ,  one_over_sqrt_3, 1.00 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }


}; // Class QuadrilateralGaussLegendreIntegrationPoints2

class QuadrilateralGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLegendreIntegrationPoints3);
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
        // Auxiliary variables for repeated terms
        const double sqrt_3_5 = std::sqrt(3.0 / 5.0);
        const double weight1 = 25.0 / 81.0;
        const double weight2 = 40.0 / 81.0;
        const double weight3 = 64.0 / 81.0;

        // Integration points
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-sqrt_3_5, -sqrt_3_5, weight1),
            IntegrationPointType( 0.0, -sqrt_3_5, weight2),
            IntegrationPointType( sqrt_3_5, -sqrt_3_5, weight1),
            IntegrationPointType(-sqrt_3_5,  0.0, weight2),
            IntegrationPointType( 0.0,  0.0, weight3),
            IntegrationPointType( sqrt_3_5,  0.0, weight2),
            IntegrationPointType(-sqrt_3_5,  sqrt_3_5, weight1),
            IntegrationPointType( 0.0,  sqrt_3_5, weight2),
            IntegrationPointType( sqrt_3_5,  sqrt_3_5, weight1)
        }};
        return s_integration_points;
    }


    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }


}; // Class QuadrilateralGaussLegendreIntegrationPoints3

class QuadrilateralGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLegendreIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 16> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()    {  return 16; }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        // Auxiliary variables for repeated terms
        const double sqrt_6_5 = std::sqrt(6.0 / 5.0);
        const double term1 = std::sqrt((3.0 + 2.0 * sqrt_6_5) / 7.0);
        const double term2 = std::sqrt((3.0 - 2.0 * sqrt_6_5) / 7.0);
        const double sqrt_5_6 = std::sqrt(5.0 / 6.0);
        const double weight1 = (0.5 - sqrt_5_6 / 6.0);
        const double weight2 = (0.5 + sqrt_5_6 / 6.0);

        // Pre-computed weights
        const double weight1_squared = weight1 * weight1;
        const double weight1_weight2 = weight1 * weight2;
        const double weight2_squared = weight2 * weight2;

        // Integration points
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-term1, -term1, weight1_squared),
            IntegrationPointType(-term1, -term2, weight1_weight2),
            IntegrationPointType(-term1, term2, weight1_weight2),
            IntegrationPointType(-term1, term1, weight1_squared),
            IntegrationPointType(-term2, -term1, weight1_weight2),
            IntegrationPointType(-term2, -term2, weight2_squared),
            IntegrationPointType(-term2, term2, weight2_squared),
            IntegrationPointType(-term2, term1, weight1_weight2),
            IntegrationPointType(term2, -term1, weight1_weight2),
            IntegrationPointType(term2, -term2, weight2_squared),
            IntegrationPointType(term2, term2, weight2_squared),
            IntegrationPointType(term2, term1, weight1_weight2),
            IntegrationPointType(term1, -term1, weight1_squared),
            IntegrationPointType(term1, -term2, weight1_weight2),
            IntegrationPointType(term1, term2, weight1_weight2),
            IntegrationPointType(term1, term1, weight1_squared)
        }};

        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }


}; // Class QuadrilateralGaussLegendreIntegrationPoints4

class QuadrilateralGaussLegendreIntegrationPoints5 {
public:
	KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLegendreIntegrationPoints5);
	typedef std::size_t SizeType;

	static const unsigned int Dimension = 2;

	typedef IntegrationPoint<2> IntegrationPointType;

	typedef std::array<IntegrationPointType, 25> IntegrationPointsArrayType;

	typedef IntegrationPointType::PointType PointType;

	static SizeType IntegrationPointsNumber() {return 25;}

    static IntegrationPointsArrayType IntegrationPoints()
    {
        static IntegrationPointsArrayType s_integration_points;
        const double a[] = {-0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664};
        const double w[] = {0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189};

        for(unsigned int i = 0; i < 5; ++i) {
            for(unsigned int j = 0; j < 5; ++j) {
                s_integration_points[5*i + j] = IntegrationPointType( a[i], a[j], w[i] * w[j]);
            }
        }

        return s_integration_points;
	}

	std::string Info() const
	{
		std::stringstream buffer;
		buffer << "Quadrilateral Gauss-Legendre quadrature 5 ";
		return buffer.str();
	}


}; // Class QuadrilateralGaussLegendreIntegrationPoints5


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.


