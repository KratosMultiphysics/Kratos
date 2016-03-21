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

//class KRATOS_API(KRATOS_CORE) PrismGaussLegendreIntegrationPoints1
//{
//public:
//    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints1);
//    typedef std::size_t SizeType;

//    static const unsigned int Dimension = 3;

//    typedef IntegrationPoint<3> IntegrationPointType;

//    typedef boost::array<IntegrationPointType, 1> IntegrationPointsArrayType;

//    typedef IntegrationPointType::PointType PointType;

//    static SizeType IntegrationPointsNumber()
//    {
//        return 1;
//    }

//    static IntegrationPointsArrayType& IntegrationPoints()
//    {

//        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5 , 0.5 );
//        return msIntegrationPoints;
//    }

//    std::string Info() const
//    {
//        std::stringstream buffer;
//        buffer << "Prism Gauss-Legendre quadrature 1 ";
//        return buffer.str();
//    }
//protected:

//private:

//    static IntegrationPointsArrayType msIntegrationPoints;

//}; // Class PrismGaussLegendreIntegrationPoints1

class KRATOS_API(KRATOS_CORE) PrismGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints1);
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

        msIntegrationPoints[0] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.16666666666666666667 );
        msIntegrationPoints[1] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.50000000000000000000 , 0.16666666666666666667 );
        msIntegrationPoints[2] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.16666666666666666667 );
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
        msIntegrationPoints[0] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.2113248654051871177 , 0.08333333333333333333 );
        msIntegrationPoints[1] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.2113248654051871177 , 0.08333333333333333333 );
        msIntegrationPoints[2] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.2113248654051871177 , 0.08333333333333333333 );
        msIntegrationPoints[3] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.7886751345948128823 , 0.08333333333333333333 );
        msIntegrationPoints[4] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.7886751345948128823 , 0.08333333333333333333 );
        msIntegrationPoints[5] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.7886751345948128823 , 0.08333333333333333333 );

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
        msIntegrationPoints[0] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.11270166537925831148 , 0.04629629629629629630 );
        msIntegrationPoints[1] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.11270166537925831148 , 0.04629629629629629630 );
        msIntegrationPoints[2] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.11270166537925831148 , 0.04629629629629629630 );
        msIntegrationPoints[3] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.07407407407407407408 );
        msIntegrationPoints[4] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.50000000000000000000 , 0.07407407407407407408 );
        msIntegrationPoints[5] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.07407407407407407408 );
        msIntegrationPoints[6] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.88729833462074168852 , 0.04629629629629629630 );
        msIntegrationPoints[7] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.88729833462074168852 , 0.04629629629629629630 );
        msIntegrationPoints[8] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.88729833462074168852 , 0.04629629629629629630 );

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

class KRATOS_API(KRATOS_CORE) PrismGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 12> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 12;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.06943184420297371239 , 0.02898790376145448812 );
        msIntegrationPoints[1]  = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.06943184420297371239 , 0.02898790376145448812 );
        msIntegrationPoints[2]  = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.06943184420297371239 , 0.02898790376145448812 );
        msIntegrationPoints[3]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.33000947820757186760 , 0.05434542957187884522 );
        msIntegrationPoints[4]  = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.33000947820757186760 , 0.05434542957187884522 );
        msIntegrationPoints[5]  = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.33000947820757186760 , 0.05434542957187884522 );
        msIntegrationPoints[6]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.66999052179242813240 , 0.05434542957187884522 );
        msIntegrationPoints[7]  = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.66999052179242813240 , 0.05434542957187884522 );
        msIntegrationPoints[8]  = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.66999052179242813240 , 0.05434542957187884522 );
        msIntegrationPoints[9]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.93056815579702628761 , 0.02898790376145448812 );
        msIntegrationPoints[10] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.93056815579702628761 , 0.02898790376145448812 );
        msIntegrationPoints[11] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.93056815579702628761 , 0.02898790376145448812 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPoints4

class KRATOS_API(KRATOS_CORE) PrismGaussLegendreIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 15> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 15;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.04691007703066800360 , 0.01974390708801575729 );
        msIntegrationPoints[1]  = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.04691007703066800360 , 0.01974390708801575729 );
        msIntegrationPoints[2]  = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.04691007703066800360 , 0.01974390708801575729 );
        msIntegrationPoints[3]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.23076534494715845448 , 0.03988572254161387234 );
        msIntegrationPoints[4]  = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.23076534494715845448 , 0.03988572254161387234 );
        msIntegrationPoints[5]  = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.23076534494715845448 , 0.03988572254161387234 );
        msIntegrationPoints[6]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.04740740740740740741 );
        msIntegrationPoints[7]  = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.50000000000000000000 , 0.04740740740740740741 );
        msIntegrationPoints[8]  = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.04740740740740740741 );
        msIntegrationPoints[9]  = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.76923465505284154552 , 0.03988572254161387234 );
        msIntegrationPoints[10] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.76923465505284154552 , 0.03988572254161387234 );
        msIntegrationPoints[11] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.76923465505284154552 , 0.03988572254161387234 );
        msIntegrationPoints[12] = IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.95308992296933199640 , 0.01974390708801575729 );
        msIntegrationPoints[13] = IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.95308992296933199640 , 0.01974390708801575729 );
        msIntegrationPoints[14] = IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.95308992296933199640 , 0.01974390708801575729 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 5 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLegendreIntegrationPoints5



///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined 


