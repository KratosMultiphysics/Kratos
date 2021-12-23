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
//                   Vicente Mataix Ferrandiz
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

class PrismGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.16666666666666666667 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.50000000000000000000 , 0.16666666666666666667 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.16666666666666666667 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPoints1

class PrismGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.2113248654051871177 , 0.08333333333333333333 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.2113248654051871177 , 0.08333333333333333333 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.2113248654051871177 , 0.08333333333333333333 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.7886751345948128823 , 0.08333333333333333333 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.7886751345948128823 , 0.08333333333333333333 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.7886751345948128823 , 0.08333333333333333333 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPoints2

class PrismGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 9> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 9;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.11270166537925831148 , 0.04629629629629629630 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.11270166537925831148 , 0.04629629629629629630 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.11270166537925831148 , 0.04629629629629629630 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.07407407407407407408 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.50000000000000000000 , 0.07407407407407407408 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.07407407407407407408 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.88729833462074168852 , 0.04629629629629629630 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.88729833462074168852 , 0.04629629629629629630 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.88729833462074168852 , 0.04629629629629629630 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPoints3

class PrismGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 12> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 12;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.06943184420297371239 , 0.02898790376145448812 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.06943184420297371239 , 0.02898790376145448812 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.06943184420297371239 , 0.02898790376145448812 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.33000947820757186760 , 0.05434542957187884522 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.33000947820757186760 , 0.05434542957187884522 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.33000947820757186760 , 0.05434542957187884522 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.66999052179242813240 , 0.05434542957187884522 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.66999052179242813240 , 0.05434542957187884522 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.66999052179242813240 , 0.05434542957187884522 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.93056815579702628761 , 0.02898790376145448812 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.93056815579702628761 , 0.02898790376145448812 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.93056815579702628761 , 0.02898790376145448812 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPoints4

class PrismGaussLegendreIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 15> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 15;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.04691007703066800360 , 0.01974390708801575729 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.04691007703066800360 , 0.01974390708801575729 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.04691007703066800360 , 0.01974390708801575729 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.23076534494715845448 , 0.03988572254161387234 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.23076534494715845448 , 0.03988572254161387234 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.23076534494715845448 , 0.03988572254161387234 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.04740740740740740741 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.50000000000000000000 , 0.04740740740740740741 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.50000000000000000000 , 0.04740740740740740741 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.76923465505284154552 , 0.03988572254161387234 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.76923465505284154552 , 0.03988572254161387234 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.76923465505284154552 , 0.03988572254161387234 ),
            IntegrationPointType( 0.66666666666666666667 , 0.16666666666666666667 , 0.95308992296933199640 , 0.01974390708801575729 ),
            IntegrationPointType( 0.16666666666666666667 , 0.66666666666666666667 , 0.95308992296933199640 , 0.01974390708801575729 ),
            IntegrationPointType( 0.16666666666666666667 , 0.16666666666666666667 , 0.95308992296933199640 , 0.01974390708801575729 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 5 ";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPoints5

/* For the formulation SPRISM we just consider one integration point in the plane */

class PrismGaussLegendreIntegrationPointsExt1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsExt1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 2> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 2;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2113248654051871177454 , 0.25000000000000000000000 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7886751345948128822546 , 0.25000000000000000000000 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 2, 1 point in plane";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPointsExt1

class PrismGaussLegendreIntegrationPointsExt2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsExt2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.1127016653792583114821 , 0.1388888888888888889 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.2222222222222222222 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.8872983346207416885180 , 0.1388888888888888889 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 3, 1 point in plane";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPointsExt2

class PrismGaussLegendreIntegrationPointsExt3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsExt3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 5> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 5;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0469100770306680036012 , 0.0592317212640472718 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2307653449471584544819 , 0.1196571676248416170 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.1422222222222222222 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7692346550528415455182 , 0.1196571676248416170 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9530899229693319963988 , 0.0592317212640472718 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 5, 1 point in plane";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPointsExt3

class PrismGaussLegendreIntegrationPointsExt4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsExt4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 7> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 7;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0254460438286207377369 , 0.0261224489795918367347 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.1292344072003027800681 , 0.069926347872319166975  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2970774243113014165467 , 0.09545751262627973624   ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.1044897959183673469388 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7029225756886985834533 , 0.09545751262627973624   ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.8707655927996972199320 , 0.069926347872319166975  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9745539561713792622631 , 0.0261224489795918367347 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 4, 1 point in plane";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPointsExt4

class PrismGaussLegendreIntegrationPointsExt5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLegendreIntegrationPointsExt5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 11> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 11;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0108856709269715035981 , 0.0139171417790434166207 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.0564687001159523504624 , 0.031395092366226156159  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.1349239972129753379533 , 0.0465725527319335628565 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.2404519353965940920372 , 0.058298441147997619980  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.3652284220238275138343 , 0.065701136127561665545  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.5000000000000000000000 , 0.0682312716944751576786 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.6347715779761724861657 , 0.065701136127561665545  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.7595480646034059079628 , 0.058298441147997619980  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.8650760027870246620467 , 0.0465725527319335628565 ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9435312998840476495376 , 0.031395092366226156159  ),
            IntegrationPointType( 0.33333333333333333333  , 0.33333333333333333333 , 0.9891143290730284964020 , 0.0139171417790434166207 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Legendre quadrature 4, 1 point in plane";
        return buffer.str();
    }


}; // Class PrismGaussLegendreIntegrationPointsExt5

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_PRISM_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined


