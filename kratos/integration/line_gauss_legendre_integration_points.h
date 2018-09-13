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

#if !defined(KRATOS_LINE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_LINE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

namespace Kratos
{
class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints1);
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
        // This is added to solve the problem of static initialization. Pooyan.
        msIntegrationPoints[0] = IntegrationPointType(0.00, 2.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints1


class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints2);
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
        msIntegrationPoints[0] = IntegrationPointType(-std::sqrt(1.00 / 3.00), 1.00);
        msIntegrationPoints[1] = IntegrationPointType( std::sqrt(1.00 / 3.00), 1.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints2


class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints3);
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
        msIntegrationPoints[0] = IntegrationPointType(-std::sqrt(3.00 / 5.00), 5.00 / 9.00);
        msIntegrationPoints[1] = IntegrationPointType( 0.00                  , 8.00 / 9.00);
        msIntegrationPoints[2] = IntegrationPointType( std::sqrt(3.00 / 5.00), 5.00 / 9.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints3



class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints4);
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
        msIntegrationPoints[0] = IntegrationPointType(-0.861136311594053, 0.347854845137454);
        msIntegrationPoints[1] = IntegrationPointType(-0.339981043584856, 0.652145154862546);
        msIntegrationPoints[2] = IntegrationPointType( 0.339981043584856, 0.652145154862546);
        msIntegrationPoints[3] = IntegrationPointType( 0.861136311594053, 0.347854845137454);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints4



class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints5);
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
        msIntegrationPoints[0] = IntegrationPointType(-0.906179845938664, 0.236926885056189);
        msIntegrationPoints[1] = IntegrationPointType(-0.538469310105683, 0.478628670499366);
        msIntegrationPoints[2] = IntegrationPointType( 0.000000000000000, 0.568888888888889);
        msIntegrationPoints[3] = IntegrationPointType( 0.538469310105683, 0.478628670499366);
        msIntegrationPoints[4] = IntegrationPointType( 0.906179845938664, 0.236926885056189);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 5 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints4



class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints6
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints6);
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
        msIntegrationPoints[0] = IntegrationPointType(-0.9324695142031521, 0.1713244923791704);
        msIntegrationPoints[1] = IntegrationPointType(-0.6612093864662645, 0.3607615730481386);
        msIntegrationPoints[2] = IntegrationPointType(-0.2386191860831969, 0.4679139345726910);
        msIntegrationPoints[3] = IntegrationPointType( 0.2386191860831969, 0.4679139345726910);
        msIntegrationPoints[4] = IntegrationPointType( 0.6612093864662645, 0.3607615730481386);
        msIntegrationPoints[5] = IntegrationPointType( 0.9324695142031521, 0.1713244923791704);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 6 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints6



class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints7
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints7);
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
        msIntegrationPoints[0] = IntegrationPointType(-0.9491079123427585, 0.1294849661688697);
        msIntegrationPoints[1] = IntegrationPointType(-0.7415311855993945, 0.2797053914892766);
        msIntegrationPoints[2] = IntegrationPointType(-0.4058451513773972, 0.3818300505051189);
        msIntegrationPoints[3] = IntegrationPointType( 0.0000000000000000, 0.4179591836734694);
        msIntegrationPoints[4] = IntegrationPointType( 0.4058451513773972, 0.3818300505051189);
        msIntegrationPoints[5] = IntegrationPointType( 0.7415311855993945, 0.2797053914892766);
        msIntegrationPoints[6] = IntegrationPointType( 0.9491079123427585, 0.1294849661688697);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 7 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints7



class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints8
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints8);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 8> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 8;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType(-0.9602898564975363, 0.1012285362903763);
        msIntegrationPoints[1] = IntegrationPointType(-0.7966664774136267, 0.2223810344533745);
        msIntegrationPoints[2] = IntegrationPointType(-0.5255324099163290, 0.3137066458778873);
        msIntegrationPoints[3] = IntegrationPointType(-0.1834346424956498, 0.3626837833783620);
        msIntegrationPoints[4] = IntegrationPointType( 0.1834346424956498, 0.3626837833783620);
        msIntegrationPoints[5] = IntegrationPointType( 0.5255324099163290, 0.3137066458778873);
        msIntegrationPoints[6] = IntegrationPointType( 0.7966664774136267, 0.2223810344533745);
        msIntegrationPoints[7] = IntegrationPointType( 0.9602898564975363, 0.1012285362903763);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 8 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints8



class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints9
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints9);
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
        msIntegrationPoints[0] = IntegrationPointType(-0.9681602395076261, 0.0812743883615744);
        msIntegrationPoints[1] = IntegrationPointType(-0.8360311073266358, 0.1806481606948574);
        msIntegrationPoints[2] = IntegrationPointType(-0.6133714327005904, 0.2606106964029354);
        msIntegrationPoints[3] = IntegrationPointType(-0.3242534234038089, 0.3123470770400029);
        msIntegrationPoints[4] = IntegrationPointType( 0.0000000000000000, 0.3302393550012598);
        msIntegrationPoints[5] = IntegrationPointType( 0.3242534234038089, 0.3123470770400029);
        msIntegrationPoints[6] = IntegrationPointType( 0.6133714327005904, 0.2606106964029354);
        msIntegrationPoints[7] = IntegrationPointType( 0.8360311073266358, 0.1806481606948574);
        msIntegrationPoints[8] = IntegrationPointType( 0.9681602395076261, 0.0812743883615744);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 9 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints9



class KRATOS_API(KRATOS_CORE) LineGaussLegendreIntegrationPoints10
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLegendreIntegrationPoints10);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 10> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 10;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType(-0.9739065285171717, 0.0666713443086881);
        msIntegrationPoints[1] = IntegrationPointType(-0.8650633666889845, 0.1494513491505806);
        msIntegrationPoints[2] = IntegrationPointType(-0.6794095682990244, 0.2190863625159820);
        msIntegrationPoints[3] = IntegrationPointType(-0.4333953941292472, 0.2692667193099963);
        msIntegrationPoints[4] = IntegrationPointType(-0.1488743389816312, 0.2955242247147529);
        msIntegrationPoints[5] = IntegrationPointType( 0.1488743389816312, 0.2955242247147529);
        msIntegrationPoints[6] = IntegrationPointType( 0.4333953941292472, 0.2692667193099963);
        msIntegrationPoints[7] = IntegrationPointType( 0.6794095682990244, 0.2190863625159820);
        msIntegrationPoints[8] = IntegrationPointType( 0.8650633666889845, 0.1494513491505806);
        msIntegrationPoints[9] = IntegrationPointType( 0.9739065285171717, 0.0666713443086881);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Legendre quadrature 10 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLegendreIntegrationPoints10


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

#endif // KRATOS_LINE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined


