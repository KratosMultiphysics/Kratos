// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
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

}; // Class PrismGaussLegendreIntegrationPointsSimple1

class PrismGaussLegendreIntegrationPointsSimple2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple2);
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

}; // Class PrismGaussLegendreIntegrationPointsSimple2

class PrismGaussLegendreIntegrationPointsSimple3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple3);
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

}; // Class PrismGaussLegendreIntegrationPointsSimple3

class PrismGaussLegendreIntegrationPointsSimple4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsSimple4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 7> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 7;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0254460438286207377369 , 0.0261224489795918367347 );
        msIntegrationPoints[1] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.1292344072003027800681 , 0.069926347872319166975  );
        msIntegrationPoints[2] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2970774243113014165467 , 0.09545751262627973624   );
        msIntegrationPoints[3] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.1044897959183673469388 );
        msIntegrationPoints[4] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7029225756886985834533 , 0.09545751262627973624   );
        msIntegrationPoints[5] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.8707655927996972199320 , 0.069926347872319166975  );
        msIntegrationPoints[6] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9745539561713792622631 , 0.0261224489795918367347 );
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

    typedef boost::array<IntegrationPointType, 11> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 11;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[ 0] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0108856709269715035981 , 0.0139171417790434166207 );
        msIntegrationPoints[ 1] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0564687001159523504624 , 0.031395092366226156159  );
        msIntegrationPoints[ 2] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.1349239972129753379533 , 0.0465725527319335628565 );
        msIntegrationPoints[ 3] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2404519353965940920372 , 0.058298441147997619980  );
        msIntegrationPoints[ 4] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.3652284220238275138343 , 0.065701136127561665545  );
        msIntegrationPoints[ 5] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.0682312716944751576786 );
        msIntegrationPoints[ 6] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.6347715779761724861657 , 0.065701136127561665545  );
        msIntegrationPoints[ 7] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7595480646034059079628 , 0.058298441147997619980  );
        msIntegrationPoints[ 8] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.8650760027870246620467 , 0.0465725527319335628565 );
        msIntegrationPoints[ 9] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9435312998840476495376 , 0.031395092366226156159  );
        msIntegrationPoints[10] = IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9891143290730284964020 , 0.0139171417790434166207 );
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

}; // Class PrismGaussLegendreIntegrationPointsSimple5

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_SIMPLE_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined
