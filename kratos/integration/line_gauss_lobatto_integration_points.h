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


#if !defined(KRATOS_LINE_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_LINE_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{
class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 1> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType(0.00, 2.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

};


class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 2> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 2;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType(-1.00, 1.00);
        msIntegrationPoints[1] = IntegrationPointType( 1.00, 1.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

};


class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints3);
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
        msIntegrationPoints[0] = IntegrationPointType(-1.00, 1.00 / 3.00);
        msIntegrationPoints[1] = IntegrationPointType( 0.00, 4.00 / 3.00);
        msIntegrationPoints[2] = IntegrationPointType( 1.00, 1.00 / 3.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints3



class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType(-1.00, 1.00 / 6.00);
        msIntegrationPoints[1] = IntegrationPointType(-std::sqrt(5.00) / 5.00, 5.00 / 6.00);
        msIntegrationPoints[2] = IntegrationPointType( std::sqrt(5.00) / 5.00, 5.00 / 6.00);
        msIntegrationPoints[3] = IntegrationPointType( 1.00, 1.00 / 6.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints4



class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints5);
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
        msIntegrationPoints[0] = IntegrationPointType(-1.00, 0.10);
        msIntegrationPoints[1] = IntegrationPointType(-std::sqrt(21.00) / 7.00, 49.00 / 90.00);
        msIntegrationPoints[2] = IntegrationPointType( 0.00, 32.00 / 45.00);
        msIntegrationPoints[3] = IntegrationPointType( std::sqrt(21.00) / 7.00, 49.00 / 90.00);
        msIntegrationPoints[4] = IntegrationPointType( 1.00, 0.10);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 5 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints4



class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints6
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints6);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType(-1.00, 1.00 / 15.00);
        msIntegrationPoints[1] = IntegrationPointType(-std::sqrt((7.00+2.00*std::sqrt(7)) / 21.00), (14.00-std::sqrt(7)) / 30.00);
        msIntegrationPoints[2] = IntegrationPointType(-std::sqrt((7.00-2.00*std::sqrt(7)) / 21.00), (14.00+std::sqrt(7)) / 30.00);
        msIntegrationPoints[3] = IntegrationPointType( std::sqrt((7.00-2.00*std::sqrt(7)) / 21.00), (14.00+std::sqrt(7)) / 30.00);
        msIntegrationPoints[4] = IntegrationPointType( std::sqrt((7.00+2.00*std::sqrt(7)) / 21.00), (14.00-std::sqrt(7)) / 30.00);
        msIntegrationPoints[5] = IntegrationPointType( 1.00, 1.00 / 15.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 6 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints6


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_LINE_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED  defined


