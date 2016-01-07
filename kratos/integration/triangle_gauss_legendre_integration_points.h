//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   JMCarbonell $
//   Date:                $Date:   December 2015 $
//   Revision:            $Revision:         1.5 $
//
//


#if !defined(KRATOS_TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) TriangleGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussLegendreIntegrationPoints1


class KRATOS_API(KRATOS_CORE) TriangleGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 );
        msIntegrationPoints[1] = IntegrationPointType( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 );
        msIntegrationPoints[2] = IntegrationPointType( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussLegendreIntegrationPoints2


class KRATOS_API(KRATOS_CORE) TriangleGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 1.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 );
        msIntegrationPoints[1] = IntegrationPointType( 3.00 / 5.00 , 1.00 / 5.00 , 25.00 / 96.00 );
        msIntegrationPoints[2] = IntegrationPointType( 1.00 / 5.00 , 3.00 / 5.00 , 25.00 / 96.00 );
        msIntegrationPoints[3] = IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -27.00 / 96.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussLegendreIntegrationPoints2


class KRATOS_API(KRATOS_CORE) TriangleGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussLegendreIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()    {  return 6; }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        const double wa = 0.054975871827661;
        const double wb = 0.1116907948390055;
        const double Na1 = 0.816847572980459;
        const double Nb1 = 0.108103018168070;
        const double Na2 = 0.091576213509771;
        const double Nb2 = 0.445948490915965;

	msIntegrationPoints[0] = IntegrationPointType( Na2, Na2, wa );
        msIntegrationPoints[1] = IntegrationPointType( Na1, Na2, wa );
        msIntegrationPoints[2] = IntegrationPointType( Na2, Na1, wa );
        msIntegrationPoints[3] = IntegrationPointType( Nb2, Nb2, wb );
        msIntegrationPoints[4] = IntegrationPointType( Nb1, Nb2, wb );
        msIntegrationPoints[5] = IntegrationPointType( Nb2, Nb1, wb );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussLegendreIntegrationPoints4


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined 


