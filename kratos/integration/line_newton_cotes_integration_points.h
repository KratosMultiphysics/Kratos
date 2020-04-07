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

#if !defined(KRATOS_LINE_NEWTON_COTES_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_LINE_NEWTON_COTES_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{
class KRATOS_API(KRATOS_CORE) LineNewtonCotesIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineNewtonCotesIntegrationPoints1);
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
        msIntegrationPoints[0]=IntegrationPointType(-0.5,1.333333333333333);
        msIntegrationPoints[1]=IntegrationPointType(0,-0.6666666666666665);
        msIntegrationPoints[2]=IntegrationPointType(0.5,1.333333333333333);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Newton-Cotes quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineNewtonCotesIntegrationPoints1


class KRATOS_API(KRATOS_CORE) LineNewtonCotesIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineNewtonCotesIntegrationPoints2);
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
        msIntegrationPoints[0]=IntegrationPointType(-0.6666666666666666,1.1);
        msIntegrationPoints[1]=IntegrationPointType(-0.3333333333333333,-1.4);
        msIntegrationPoints[2]=IntegrationPointType(0,2.6);
        msIntegrationPoints[3]=IntegrationPointType(0.3333333333333333,-1.4);
        msIntegrationPoints[4]=IntegrationPointType(0.6666666666666666,1.1);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Newton-Cotes quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineNewtonCotesIntegrationPoints2


class KRATOS_API(KRATOS_CORE) LineNewtonCotesIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineNewtonCotesIntegrationPoints3);
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
        msIntegrationPoints[0]=IntegrationPointType(-0.75,0.9735449735449742);
        msIntegrationPoints[1]=IntegrationPointType(-0.5,-2.019047619047615);
        msIntegrationPoints[2]=IntegrationPointType(-0.25,4.647619047619042);
        msIntegrationPoints[3]=IntegrationPointType(0,-5.204232804232804);
        msIntegrationPoints[4]=IntegrationPointType(0.25,4.647619047619049);
        msIntegrationPoints[5]=IntegrationPointType(0.5,-2.019047619047616);
        msIntegrationPoints[6]=IntegrationPointType(0.75,0.9735449735449739);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Newton-Cotes quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineNewtonCotesIntegrationPoints3



class KRATOS_API(KRATOS_CORE) LineNewtonCotesIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineNewtonCotesIntegrationPoints4);
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
        msIntegrationPoints[0]=IntegrationPointType(-0.8,0.8917548500881828);
        msIntegrationPoints[1]=IntegrationPointType(-0.6,-2.577160493827184);
        msIntegrationPoints[2]=IntegrationPointType(-0.4,7.350088183421553);
        msIntegrationPoints[3]=IntegrationPointType(-0.2,-12.14065255731907);
        msIntegrationPoints[4]=IntegrationPointType(0,14.95194003527322);
        msIntegrationPoints[5]=IntegrationPointType(0.2,-12.14065255731914);
        msIntegrationPoints[6]=IntegrationPointType(0.4,7.350088183421514);
        msIntegrationPoints[7]=IntegrationPointType(0.6,-2.577160493827156);
        msIntegrationPoints[8]=IntegrationPointType(0.8,0.8917548500881831);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Newton-Cotes quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineNewtonCotesIntegrationPoints4



class KRATOS_API(KRATOS_CORE) LineNewtonCotesIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineNewtonCotesIntegrationPoints5);
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
        msIntegrationPoints[0]=IntegrationPointType(-0.8333333333333334,0.8334199134199013);
        msIntegrationPoints[1]=IntegrationPointType(-0.6666666666666666,-3.097056277056003);
        msIntegrationPoints[2]=IntegrationPointType(-0.5,10.65437229437228);
        msIntegrationPoints[3]=IntegrationPointType(-0.3333333333333333,-23.0561038961028);
        msIntegrationPoints[4]=IntegrationPointType(-0.1666666666666667,37.05246753246604);
        msIntegrationPoints[5]=IntegrationPointType(0,-42.77419913420147);
        msIntegrationPoints[6]=IntegrationPointType(0.1666666666666667,37.05246753246799);
        msIntegrationPoints[7]=IntegrationPointType(0.3333333333333333,-23.05610389610366);
        msIntegrationPoints[8]=IntegrationPointType(0.5,10.65437229437222);
        msIntegrationPoints[9]=IntegrationPointType(0.6666666666666666,-3.097056277056263);
        msIntegrationPoints[10]=IntegrationPointType(0.8333333333333334,0.8334199134199121);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Newton-Cotes quadrature 5 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineNewtonCotesIntegrationPoints5


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

#endif // KRATOS_LINE_NEWTON_COTES_INTEGRATION_POINTS_H_INCLUDED  defined


