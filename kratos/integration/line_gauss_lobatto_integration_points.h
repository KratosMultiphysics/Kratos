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

}; // Class LineGaussLobattoIntegrationPoints5



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



class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints7
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints7);
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
        msIntegrationPoints[0] = IntegrationPointType(-1.00, 1.00 / 21.00);
        msIntegrationPoints[1] = IntegrationPointType(-std::sqrt((5.00/11.00) + (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 - 7.00*std::sqrt(15.00)) / 350.00);
        msIntegrationPoints[2] = IntegrationPointType(-std::sqrt((5.00/11.00) - (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 + 7.00*std::sqrt(15.00)) / 350.00);
        msIntegrationPoints[3] = IntegrationPointType(0.00, 256.00/525.00);
        msIntegrationPoints[4] = IntegrationPointType(std::sqrt((5.00/11.00) - (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 + 7.00*std::sqrt(15.00)) / 350.00);
        msIntegrationPoints[5] = IntegrationPointType(std::sqrt((5.00/11.00) + (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 - 7.00*std::sqrt(15.00)) / 350.00);
        msIntegrationPoints[6] = IntegrationPointType(-1.00, 1.00 / 21.00);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 7 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints7



class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints8
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints8);
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
        msIntegrationPoints[0] = IntegrationPointType(-1.00             , 0.035714285714286);
        msIntegrationPoints[1] = IntegrationPointType(-0.871740148509607, 0.210704227143506);
        msIntegrationPoints[2] = IntegrationPointType(-0.591700181433142, 0.341122692483504);
        msIntegrationPoints[3] = IntegrationPointType(-0.209299217902479, 0.412458794658704);
        msIntegrationPoints[4] = IntegrationPointType( 0.209299217902479, 0.412458794658704);
        msIntegrationPoints[5] = IntegrationPointType( 0.591700181433142, 0.341122692483504);
        msIntegrationPoints[6] = IntegrationPointType( 0.871740148509607, 0.210704227143506);
        msIntegrationPoints[7] = IntegrationPointType( 1.00             , 0.035714285714286);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 8 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints8



class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints9
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints9);
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
        msIntegrationPoints[0] = IntegrationPointType(-1.00             , 0.027777777777778);
        msIntegrationPoints[1] = IntegrationPointType(-0.899757995411460, 0.165495361560806);
        msIntegrationPoints[2] = IntegrationPointType(-0.677186279510738, 0.274538712500162);
        msIntegrationPoints[3] = IntegrationPointType(-0.363117463826178, 0.346428510973046);
        msIntegrationPoints[4] = IntegrationPointType( 0.00             , 0.371519274376417);
        msIntegrationPoints[5] = IntegrationPointType( 0.363117463826178, 0.346428510973046);
        msIntegrationPoints[6] = IntegrationPointType( 0.677186279510738, 0.274538712500162);
        msIntegrationPoints[7] = IntegrationPointType( 0.899757995411460, 0.165495361560806);
        msIntegrationPoints[8] = IntegrationPointType( 1.00             , 0.027777777777778);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 9 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints9



class KRATOS_API(KRATOS_CORE) LineGaussLobattoIntegrationPoints10
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LineGaussLobattoIntegrationPoints10);
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
        msIntegrationPoints[0] = IntegrationPointType(-1.00             , 0.022222222222222);
        msIntegrationPoints[1] = IntegrationPointType(-0.919533908166459, 0.133305990851070);
        msIntegrationPoints[2] = IntegrationPointType(-0.738773865105505, 0.224889342063126);
        msIntegrationPoints[3] = IntegrationPointType(-0.477924949810444, 0.292042683679684);
        msIntegrationPoints[4] = IntegrationPointType(-0.165278957666387, 0.327539761183897);
        msIntegrationPoints[5] = IntegrationPointType( 0.165278957666387, 0.327539761183897);
        msIntegrationPoints[6] = IntegrationPointType( 0.477924949810444, 0.292042683679684);
        msIntegrationPoints[7] = IntegrationPointType( 0.738773865105505, 0.224889342063126);
        msIntegrationPoints[8] = IntegrationPointType( 0.919533908166459, 0.133305990851070);
        msIntegrationPoints[9] = IntegrationPointType( 1.00             , 0.022222222222222);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 10 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class LineGaussLobattoIntegrationPoints10


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_LINE_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED  defined


