//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_QUADRILATERAL_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{

class QuadrilateralCollocationIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints1);
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
            IntegrationPointType(-0.500000000000,-0.500000000000,1.00000000000),
            IntegrationPointType(-0.500000000000,0.500000000000,1.00000000000),
            IntegrationPointType(0.500000000000,-0.500000000000,1.00000000000),
            IntegrationPointType(0.500000000000,0.500000000000,1.00000000000)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 1 ";
        return buffer.str();
    }


}; // Class QuadrilateralCollocationIntegrationPoints1

class QuadrilateralCollocationIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints2);
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
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.666666666667,-0.666666666667,0.444444444444),
            IntegrationPointType(-0.666666666667,0,0.444444444444),
            IntegrationPointType(-0.666666666667,0.666666666667,0.444444444444),
            IntegrationPointType(0,-0.666666666667,0.444444444444),
            IntegrationPointType(0,0,0.444444444444),
            IntegrationPointType(0,0.666666666667,0.444444444444),
            IntegrationPointType(0.666666666667,-0.666666666667,0.444444444444),
            IntegrationPointType(0.666666666667,0,0.444444444444),
            IntegrationPointType(0.666666666667,0.666666666667,0.444444444444)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 2 ";
        return buffer.str();
    }


}; // Class QuadrilateralCollocationIntegrationPoints2

class QuadrilateralCollocationIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 16> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 16;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.750000000000,-0.750000000000,0.250000000000),
            IntegrationPointType(-0.750000000000,-0.250000000000,0.250000000000),
            IntegrationPointType(-0.750000000000,0.250000000000,0.250000000000),
            IntegrationPointType(-0.750000000000,0.750000000000,0.250000000000),
            IntegrationPointType(-0.250000000000,-0.750000000000,0.250000000000),
            IntegrationPointType(-0.250000000000,-0.250000000000,0.250000000000),
            IntegrationPointType(-0.250000000000,0.250000000000,0.250000000000),
            IntegrationPointType(-0.250000000000,0.750000000000,0.250000000000),
            IntegrationPointType(0.250000000000,-0.750000000000,0.250000000000),
            IntegrationPointType(0.250000000000,-0.250000000000,0.250000000000),
            IntegrationPointType(0.250000000000,0.250000000000,0.250000000000),
            IntegrationPointType(0.250000000000,0.750000000000,0.250000000000),
            IntegrationPointType(0.750000000000,-0.750000000000,0.250000000000),
            IntegrationPointType(0.750000000000,-0.250000000000,0.250000000000),
            IntegrationPointType(0.750000000000,0.250000000000,0.250000000000),
            IntegrationPointType(0.750000000000,0.750000000000,0.250000000000)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 3 ";
        return buffer.str();
    }


}; // Class QuadrilateralCollocationIntegrationPoints3

class QuadrilateralCollocationIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 25> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 25;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.800000000000,-0.800000000000,0.160000000000),
            IntegrationPointType(-0.800000000000,-0.400000000000,0.160000000000),
            IntegrationPointType(-0.800000000000, 0.000000000000,0.160000000000),
            IntegrationPointType(-0.800000000000,0.400000000000,0.160000000000),
            IntegrationPointType(-0.800000000000,0.800000000000,0.160000000000),
            IntegrationPointType(-0.400000000000,-0.800000000000,0.160000000000),
            IntegrationPointType(-0.400000000000,-0.400000000000,0.160000000000),
            IntegrationPointType(-0.400000000000,0.000000000000,0.160000000000),
            IntegrationPointType(-0.400000000000,0.400000000000,0.160000000000),
            IntegrationPointType(-0.400000000000,0.800000000000,0.160000000000),
            IntegrationPointType(0.0000000000000,-0.800000000000,0.160000000000),
            IntegrationPointType(0.0000000000000,-0.400000000000,0.160000000000),
            IntegrationPointType(0.0000000000000, 0.000000000000,0.160000000000),
            IntegrationPointType(0.0000000000000, 0.400000000000,0.160000000000),
            IntegrationPointType(0.0000000000000,0.800000000000,0.160000000000),
            IntegrationPointType(0.400000000000,-0.800000000000,0.160000000000),
            IntegrationPointType(0.400000000000,-0.400000000000,0.160000000000),
            IntegrationPointType(0.400000000000, 0.000000000000,0.160000000000),
            IntegrationPointType(0.400000000000,0.400000000000,0.160000000000),
            IntegrationPointType(0.400000000000,0.800000000000,0.160000000000),
            IntegrationPointType(0.800000000000,-0.800000000000,0.160000000000),
            IntegrationPointType(0.800000000000,-0.400000000000,0.160000000000),
            IntegrationPointType(0.800000000000, 0.000000000000,0.160000000000),
            IntegrationPointType(0.800000000000,0.400000000000,0.160000000000),
            IntegrationPointType(0.800000000000,0.800000000000,0.160000000000)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 4 ";
        return buffer.str();
    }


}; // Class QuadrilateralCollocationIntegrationPoints4

class QuadrilateralCollocationIntegrationPoints5 {
public:
	KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints5);
	typedef std::size_t SizeType;

	static const unsigned int Dimension = 2;

	typedef IntegrationPoint<2> IntegrationPointType;

	typedef std::array<IntegrationPointType, 36> IntegrationPointsArrayType;

	typedef IntegrationPointType::PointType PointType;

	static SizeType IntegrationPointsNumber()
        {
            return 36;
        }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-0.833333333333,-0.833333333333,0.111111111111),
            IntegrationPointType(-0.833333333333,-0.500000000000,0.111111111111),
            IntegrationPointType(-0.833333333333,-0.166666666667,0.111111111111),
            IntegrationPointType(-0.833333333333,0.166666666667,0.111111111111),
            IntegrationPointType(-0.833333333333,0.500000000000,0.111111111111),
            IntegrationPointType(-0.833333333333,0.833333333333,0.111111111111),
            IntegrationPointType(-0.500000000000,-0.833333333333,0.111111111111),
            IntegrationPointType(-0.500000000000,-0.500000000000,0.111111111111),
            IntegrationPointType(-0.500000000000,-0.166666666667,0.111111111111),
            IntegrationPointType(-0.500000000000,0.166666666667,0.111111111111),
            IntegrationPointType(-0.500000000000,0.500000000000,0.111111111111),
            IntegrationPointType(-0.500000000000,0.833333333333,0.111111111111),
            IntegrationPointType(-0.166666666667,-0.833333333333,0.111111111111),
            IntegrationPointType(-0.166666666667,-0.500000000000,0.111111111111),
            IntegrationPointType(-0.166666666667,-0.166666666667,0.111111111111),
            IntegrationPointType(-0.166666666667,0.166666666667,0.111111111111),
            IntegrationPointType(-0.166666666667,0.500000000000,0.111111111111),
            IntegrationPointType(-0.166666666667,0.833333333333,0.111111111111),
            IntegrationPointType(0.166666666667,-0.833333333333,0.111111111111),
            IntegrationPointType(0.166666666667,-0.500000000000,0.111111111111),
            IntegrationPointType(0.166666666667,-0.166666666667,0.111111111111),
            IntegrationPointType(0.166666666667,0.166666666667,0.111111111111),
            IntegrationPointType(0.166666666667,0.500000000000,0.111111111111),
            IntegrationPointType(0.166666666667,0.833333333333,0.111111111111),
            IntegrationPointType(0.500000000000,-0.833333333333,0.111111111111),
            IntegrationPointType(0.500000000000,-0.500000000000,0.111111111111),
            IntegrationPointType(0.500000000000,-0.166666666667,0.111111111111),
            IntegrationPointType(0.500000000000,0.166666666667,0.111111111111),
            IntegrationPointType(0.500000000000,0.500000000000,0.111111111111),
            IntegrationPointType(0.500000000000,0.833333333333,0.111111111111),
            IntegrationPointType(0.833333333333,-0.833333333333,0.111111111111),
            IntegrationPointType(0.833333333333,-0.500000000000,0.111111111111),
            IntegrationPointType(0.833333333333,-0.166666666667,0.111111111111),
            IntegrationPointType(0.833333333333,0.166666666667,0.111111111111),
            IntegrationPointType(0.833333333333,0.500000000000,0.111111111111),
            IntegrationPointType(0.833333333333,0.833333333333,0.111111111111)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 5 ";
        return buffer.str();
    }


}; // Class QuadrilateralCollocationIntegrationPoints5


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_QUADRILATERAL_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED  defined


