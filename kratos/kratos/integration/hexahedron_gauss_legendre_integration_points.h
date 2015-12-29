//
//   Project Name:        Kratos
//   Last Modified by:    $Author :    JMCarbonell $
//   Date:                $Date:     December 2015 $
//   Revision:            $Revision:           1.5 $
//
//


#if !defined(KRATOS_HEXAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_HEXAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.00 , 0.00, 0.00 , 8.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints1

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 8> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 8;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        msIntegrationPoints[1] = IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        msIntegrationPoints[2] = IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        msIntegrationPoints[3] = IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
        msIntegrationPoints[4] = IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
        msIntegrationPoints[5] = IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
        msIntegrationPoints[6] = IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
        msIntegrationPoints[7] = IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints2

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 27> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 27;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[ 0] = IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );
        msIntegrationPoints[ 1] = IntegrationPointType(                   0.0 , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 200.00/729.00 );
        msIntegrationPoints[ 2] = IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );

        msIntegrationPoints[ 3] = IntegrationPointType( -std::sqrt(3.00/5.00) ,                   0.0, -std::sqrt(3.00/5.00), 200.00/729.00 );
        msIntegrationPoints[ 4] = IntegrationPointType(                   0.0 ,                   0.0, -std::sqrt(3.00/5.00), 320.00/729.00 );
        msIntegrationPoints[ 5] = IntegrationPointType(  std::sqrt(3.00/5.00) ,                   0.0, -std::sqrt(3.00/5.00), 200.00/729.00 );

        msIntegrationPoints[ 6] = IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );
        msIntegrationPoints[ 7] = IntegrationPointType(                   0.0 ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 200.00/729.00 );
        msIntegrationPoints[ 8] = IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );

        msIntegrationPoints[ 9] = IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );
        msIntegrationPoints[10] = IntegrationPointType(                   0.0 , -std::sqrt(3.00/5.00),                   0.0, 320.00/729.00 );
        msIntegrationPoints[11] = IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );

        msIntegrationPoints[12] = IntegrationPointType( -std::sqrt(3.00/5.00) ,                   0.0,                   0.0, 320.00/729.00 );
        msIntegrationPoints[13] = IntegrationPointType(                   0.0 ,                   0.0,                   0.0, 512.00/729.00 );
        msIntegrationPoints[14] = IntegrationPointType(  std::sqrt(3.00/5.00) ,                   0.0,                   0.0, 320.00/729.00 );

        msIntegrationPoints[15] = IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );
        msIntegrationPoints[16] = IntegrationPointType(                   0.0 ,  std::sqrt(3.00/5.00),                   0.0, 320.00/729.00 );
        msIntegrationPoints[17] = IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );

        msIntegrationPoints[18] = IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );
        msIntegrationPoints[19] = IntegrationPointType(                   0.0 , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 200.00/729.00 );
        msIntegrationPoints[20] = IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );

        msIntegrationPoints[21] = IntegrationPointType( -std::sqrt(3.00/5.00) ,                   0.0,  std::sqrt(3.00/5.00), 200.00/729.00 );
        msIntegrationPoints[22] = IntegrationPointType(                   0.0 ,                   0.0,  std::sqrt(3.00/5.00), 320.00/729.00 );
        msIntegrationPoints[23] = IntegrationPointType(  std::sqrt(3.00/5.00) ,                   0.0,  std::sqrt(3.00/5.00), 200.00/729.00 );

        msIntegrationPoints[24] = IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );
        msIntegrationPoints[25] = IntegrationPointType(                   0.0 ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 200.00/729.00 );
        msIntegrationPoints[26] = IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexadra Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints3

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_HEXAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined 


