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
#if !defined(KRATOS_TRIANGLE_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_TRIANGLE_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) TriangleCollocationIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleCollocationIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(0.166666666667,0.166666666667,0.166666666667);
        msIntegrationPoints[1]=IntegrationPointType(0.166666666667,0.666666666667,0.166666666667);
        msIntegrationPoints[2]=IntegrationPointType(0.666666666667,0.666666666667,0.166666666667);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Collocation quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleCollocationIntegrationPoints1


class KRATOS_API(KRATOS_CORE) TriangleCollocationIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleCollocationIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(0.111111111111,0.111111111111,0.0833333333333);
        msIntegrationPoints[1]=IntegrationPointType(0.111111111111,0.444444444444,0.0833333333333);
        msIntegrationPoints[2]=IntegrationPointType(0.111111111111,0.777777777778,0.0833333333333);
        msIntegrationPoints[3]=IntegrationPointType(0.444444444444,0.444444444444,0.0833333333333);
        msIntegrationPoints[4]=IntegrationPointType(0.444444444444,0.777777777778,0.0833333333333);
        msIntegrationPoints[5]=IntegrationPointType(0.777777777778,0.777777777778,0.0833333333333);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Collocation quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleCollocationIntegrationPoints2


class KRATOS_API(KRATOS_CORE) TriangleCollocationIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleCollocationIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 10> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 10;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(0.0833333333333,0.0833333333333,0.0500000000000);
        msIntegrationPoints[1]=IntegrationPointType(0.0833333333333,0.333333333333,0.0500000000000);
        msIntegrationPoints[2]=IntegrationPointType(0.0833333333333,0.583333333333,0.0500000000000);
        msIntegrationPoints[3]=IntegrationPointType(0.0833333333333,0.833333333333,0.0500000000000);
        msIntegrationPoints[4]=IntegrationPointType(0.333333333333,0.333333333333,0.0500000000000);
        msIntegrationPoints[5]=IntegrationPointType(0.333333333333,0.583333333333,0.0500000000000);
        msIntegrationPoints[6]=IntegrationPointType(0.333333333333,0.833333333333,0.0500000000000);
        msIntegrationPoints[7]=IntegrationPointType(0.583333333333,0.583333333333,0.0500000000000);
        msIntegrationPoints[8]=IntegrationPointType(0.583333333333,0.833333333333,0.0500000000000);
        msIntegrationPoints[9]=IntegrationPointType(0.833333333333,0.833333333333,0.0500000000000);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Collocation quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleCollocationIntegrationPoints2


class KRATOS_API(KRATOS_CORE) TriangleCollocationIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleCollocationIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 15> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 15;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(0.0666666666667,0.0666666666667,0.0333333333333);
        msIntegrationPoints[1]=IntegrationPointType(0.0666666666667,0.266666666667,0.0333333333333);
        msIntegrationPoints[2]=IntegrationPointType(0.0666666666667,0.466666666667,0.0333333333333);
        msIntegrationPoints[3]=IntegrationPointType(0.0666666666667,0.666666666667,0.0333333333333);
        msIntegrationPoints[4]=IntegrationPointType(0.0666666666667,0.866666666667,0.0333333333333);
        msIntegrationPoints[5]=IntegrationPointType(0.266666666667,0.266666666667,0.0333333333333);
        msIntegrationPoints[6]=IntegrationPointType(0.266666666667,0.466666666667,0.0333333333333);
        msIntegrationPoints[7]=IntegrationPointType(0.266666666667,0.666666666667,0.0333333333333);
        msIntegrationPoints[8]=IntegrationPointType(0.266666666667,0.866666666667,0.0333333333333);
        msIntegrationPoints[9]=IntegrationPointType(0.466666666667,0.466666666667,0.0333333333333);
        msIntegrationPoints[10]=IntegrationPointType(0.466666666667,0.666666666667,0.0333333333333);
        msIntegrationPoints[11]=IntegrationPointType(0.466666666667,0.866666666667,0.0333333333333);
        msIntegrationPoints[12]=IntegrationPointType(0.666666666667,0.666666666667,0.0333333333333);
        msIntegrationPoints[13]=IntegrationPointType(0.666666666667,0.866666666667,0.0333333333333);
        msIntegrationPoints[14]=IntegrationPointType(0.866666666667,0.866666666667,0.0333333333333);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Collocation quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleCollocationIntegrationPoints4

class KRATOS_API(KRATOS_CORE) TriangleCollocationIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleCollocationIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 21> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 21;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(0.0555555555556,0.0555555555556,0.0238095238095);
        msIntegrationPoints[1]=IntegrationPointType(0.0555555555556,0.222222222222,0.0238095238095);
        msIntegrationPoints[2]=IntegrationPointType(0.0555555555556,0.388888888889,0.0238095238095);
        msIntegrationPoints[3]=IntegrationPointType(0.0555555555556,0.555555555556,0.0238095238095);
        msIntegrationPoints[4]=IntegrationPointType(0.0555555555556,0.722222222222,0.0238095238095);
        msIntegrationPoints[5]=IntegrationPointType(0.0555555555556,0.888888888889,0.0238095238095);
        msIntegrationPoints[6]=IntegrationPointType(0.222222222222,0.222222222222,0.0238095238095);
        msIntegrationPoints[7]=IntegrationPointType(0.222222222222,0.388888888889,0.0238095238095);
        msIntegrationPoints[8]=IntegrationPointType(0.222222222222,0.555555555556,0.0238095238095);
        msIntegrationPoints[9]=IntegrationPointType(0.222222222222,0.722222222222,0.0238095238095);
        msIntegrationPoints[10]=IntegrationPointType(0.222222222222,0.888888888889,0.0238095238095);
        msIntegrationPoints[11]=IntegrationPointType(0.388888888889,0.388888888889,0.0238095238095);
        msIntegrationPoints[12]=IntegrationPointType(0.388888888889,0.555555555556,0.0238095238095);
        msIntegrationPoints[13]=IntegrationPointType(0.388888888889,0.722222222222,0.0238095238095);
        msIntegrationPoints[14]=IntegrationPointType(0.388888888889,0.888888888889,0.0238095238095);
        msIntegrationPoints[15]=IntegrationPointType(0.555555555556,0.555555555556,0.0238095238095);
        msIntegrationPoints[16]=IntegrationPointType(0.555555555556,0.722222222222,0.0238095238095);
        msIntegrationPoints[17]=IntegrationPointType(0.555555555556,0.888888888889,0.0238095238095);
        msIntegrationPoints[18]=IntegrationPointType(0.722222222222,0.722222222222,0.0238095238095);
        msIntegrationPoints[19]=IntegrationPointType(0.722222222222,0.888888888889,0.0238095238095);
        msIntegrationPoints[20]=IntegrationPointType(0.888888888889,0.888888888889,0.0238095238095);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Collocation quadrature 5 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleCollocationIntegrationPoints5


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_TRIANGLE_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED  defined


