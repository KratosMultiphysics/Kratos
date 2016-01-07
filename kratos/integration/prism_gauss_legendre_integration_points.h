//
//   Project Name:        Kratos
//   Last Modified by:    $Author: JMCarbonell/VMataix $
//   Date:                $Date:          January 2016 $
//   Revision:            $Revision:               1.4 $
//
//


#if !defined(KRATOS_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) PrismGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints1);
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

        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.00 , 0.5 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPoints1

class KRATOS_API(KRATOS_CORE) PrismGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , -0.57735026918962576451 , 0.08333333333333333333 );
        msIntegrationPoints[1] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , -0.57735026918962576451 , 0.08333333333333333333 );
        msIntegrationPoints[2] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , -0.57735026918962576451 , 0.08333333333333333333 );
        msIntegrationPoints[3] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 ,  0.57735026918962576451 , 0.08333333333333333333 );
        msIntegrationPoints[4] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 ,  0.57735026918962576451 , 0.08333333333333333333 );
        msIntegrationPoints[5] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 ,  0.57735026918962576451 , 0.08333333333333333333 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPoints2

class KRATOS_API(KRATOS_CORE) PrismGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 9> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 9;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , -0.77459666924148337704 , 0.04629629629629629630 );
        msIntegrationPoints[1] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , -0.77459666924148337704 , 0.04629629629629629630 );
        msIntegrationPoints[2] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , -0.77459666924148337704 , 0.04629629629629629630 );
        msIntegrationPoints[3] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 ,  0.00000000000000000000 , 0.07407407407407407408 );
        msIntegrationPoints[4] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 ,  0.00000000000000000000 , 0.07407407407407407408 );
        msIntegrationPoints[5] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 ,  0.00000000000000000000 , 0.07407407407407407408 );
        msIntegrationPoints[6] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 ,  0.77459666924148337704 , 0.04629629629629629630 );
        msIntegrationPoints[7] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 ,  0.77459666924148337704 , 0.04629629629629629630 );
        msIntegrationPoints[8] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 ,  0.77459666924148337704 , 0.04629629629629629630 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPoints3

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined 


