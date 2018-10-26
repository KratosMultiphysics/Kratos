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


#if !defined(KRATOS_LINE_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_LINE_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{
class KRATOS_API(KRATOS_CORE) LineCollocationIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineCollocationIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        // This is added to solve the problem of static initialization. Pooyan.
        msIntegrationPoints[0]=IntegrationPointType(-0.666666666667,0.666666666667);
        msIntegrationPoints[1]=IntegrationPointType(0,0.666666666667);
        msIntegrationPoints[2]=IntegrationPointType(0.666666666667,0.666666666667);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Collocation quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineCollocationIntegrationPoints1


class KRATOS_API(KRATOS_CORE) LineCollocationIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineCollocationIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 5> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 5;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(-0.800000000000,0.400000000000);
        msIntegrationPoints[1]=IntegrationPointType(-0.400000000000,0.400000000000);
        msIntegrationPoints[2]=IntegrationPointType(0.000000000000,0.400000000000);
        msIntegrationPoints[3]=IntegrationPointType(0.400000000000,0.400000000000);
        msIntegrationPoints[4]=IntegrationPointType(0.800000000000,0.400000000000);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Collocation quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineCollocationIntegrationPoints2


class KRATOS_API(KRATOS_CORE) LineCollocationIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineCollocationIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 7> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 7;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(-0.857142857143,0.285714285714);
        msIntegrationPoints[1]=IntegrationPointType(-0.571428571429,0.285714285714);
        msIntegrationPoints[2]=IntegrationPointType(-0.285714285714,0.285714285714);
        msIntegrationPoints[3]=IntegrationPointType(0,0.285714285714);
        msIntegrationPoints[4]=IntegrationPointType(0.285714285714,0.285714285714);
        msIntegrationPoints[5]=IntegrationPointType(0.571428571429,0.285714285714);
        msIntegrationPoints[6]=IntegrationPointType(0.857142857143,0.285714285714);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Collocation quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineCollocationIntegrationPoints3



class KRATOS_API(KRATOS_CORE) LineCollocationIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineCollocationIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 9> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 9;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(-0.888888888889,0.222222222222);
        msIntegrationPoints[1]=IntegrationPointType(-0.666666666667,0.222222222222);
        msIntegrationPoints[2]=IntegrationPointType(-0.444444444444,0.222222222222);
        msIntegrationPoints[3]=IntegrationPointType(-0.222222222222,0.222222222222);
        msIntegrationPoints[4]=IntegrationPointType(0,0.222222222222);
        msIntegrationPoints[5]=IntegrationPointType(0.222222222222,0.222222222222);
        msIntegrationPoints[6]=IntegrationPointType(0.444444444444,0.222222222222);
        msIntegrationPoints[7]=IntegrationPointType(0.666666666667,0.222222222222);
        msIntegrationPoints[8]=IntegrationPointType(0.888888888889,0.222222222222);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Collocation quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineCollocationIntegrationPoints4



class KRATOS_API(KRATOS_CORE) LineCollocationIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineCollocationIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 11> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 11;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(-0.909090909091,0.181818181818);
        msIntegrationPoints[1]=IntegrationPointType(-0.727272727273,0.181818181818);
        msIntegrationPoints[2]=IntegrationPointType(-0.545454545455,0.181818181818);
        msIntegrationPoints[3]=IntegrationPointType(-0.363636363636,0.181818181818);
        msIntegrationPoints[4]=IntegrationPointType(-0.181818181818,0.181818181818);
        msIntegrationPoints[5]=IntegrationPointType(0,0.181818181818);
        msIntegrationPoints[6]=IntegrationPointType(0.181818181818,0.181818181818);
        msIntegrationPoints[7]=IntegrationPointType(0.363636363636,0.181818181818);
        msIntegrationPoints[8]=IntegrationPointType(0.545454545455,0.181818181818);
        msIntegrationPoints[9]=IntegrationPointType(0.727272727273,0.181818181818);
        msIntegrationPoints[10]=IntegrationPointType(0.909090909091,0.181818181818);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Collocation quadrature 5 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineCollocationIntegrationPoints5


///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_LINE_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED  defined


