//
//   Project Name:        Kratos
//   Last Modified by:    $Author: VMataix $
//   Date:                $Date:          January 2016 $
//   Revision:            $Revision:               0.0 $
//
//

#if !defined(KRATOS_SIMPLE_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_SIMPLE_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{

/* For the formulation SPRISM we just consider one integration point in the plane */

class PrismGaussLegendreIntegrationPointsSimple1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple1);
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

        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.50000000000000000000000 , 0.50000000000000000000000 );
        return msIntegrationPoints;
    }
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 1, 1 point in plane";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPointsSimple1

class PrismGaussLegendreIntegrationPointsSimple2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 2> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 2;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {

        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2113248654051871177454 , 0.25000000000000000000000 );
        msIntegrationPoints[1] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7886751345948128822546 , 0.25000000000000000000000 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 2, 1 point in plane";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPointsSimple2

class PrismGaussLegendreIntegrationPointsSimple3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.1127016653792583114821 , 0.1388888888888888889 );
        msIntegrationPoints[1] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.2222222222222222222 );
        msIntegrationPoints[2] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.8872983346207416885180 , 0.1388888888888888889 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 3, 1 point in plane";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPointsSimple3

class PrismGaussLegendreIntegrationPointsSimple4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0694318442029737123881 , 0.0869637112843634643 );
        msIntegrationPoints[1] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.3300094782075718675987 , 0.1630362887156365357 );
        msIntegrationPoints[2] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.6699905217924281324013 , 0.1630362887156365357 );
        msIntegrationPoints[3] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9305681557970262876120 , 0.0869637112843634643 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 4, 1 point in plane";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPointsSimple4

class PrismGaussLegendreIntegrationPointsSimple5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 5> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 5;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0469100770306680036012 , 0.0592317212640472718 );
        msIntegrationPoints[1] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2307653449471584544819 , 0.1196571676248416170 );
        msIntegrationPoints[2] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.1422222222222222222 );
        msIntegrationPoints[3] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7692346550528415455182 , 0.1196571676248416170 );
        msIntegrationPoints[4] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9530899229693319963988 , 0.0592317212640472718 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 5, 1 point in plane";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPointsSimple5

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_SIMPLE_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined
