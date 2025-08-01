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


#pragma once


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{
class SolidMechanicsLineGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 1> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(0.00, 2.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }


};


class SolidMechanicsLineGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 2> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 2;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00, 1.00),
            IntegrationPointType( 1.00, 1.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 2 ";
        return buffer.str();
    }


};


class SolidMechanicsLineGaussLobattoIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00, 1.00 / 3.00),
            IntegrationPointType( 0.00, 4.00 / 3.00),
            IntegrationPointType( 1.00, 1.00 / 3.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 3 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints3



class SolidMechanicsLineGaussLobattoIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00, 1.00 / 6.00),
            IntegrationPointType(-std::sqrt(5.00) / 5.00, 5.00 / 6.00),
            IntegrationPointType( std::sqrt(5.00) / 5.00, 5.00 / 6.00),
            IntegrationPointType( 1.00, 1.00 / 6.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 4 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints4



class SolidMechanicsLineGaussLobattoIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 5> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 5;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00, 0.10),
            IntegrationPointType(-std::sqrt(21.00) / 7.00, 49.00 / 90.00),
            IntegrationPointType( 0.00, 32.00 / 45.00),
            IntegrationPointType( std::sqrt(21.00) / 7.00, 49.00 / 90.00),
            IntegrationPointType( 1.00, 0.10)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 5 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints5



class SolidMechanicsLineGaussLobattoIntegrationPoints6
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints6);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00, 1.00 / 15.00),
            IntegrationPointType(-std::sqrt((7.00+2.00*std::sqrt(7)) / 21.00), (14.00-std::sqrt(7)) / 30.00),
            IntegrationPointType(-std::sqrt((7.00-2.00*std::sqrt(7)) / 21.00), (14.00+std::sqrt(7)) / 30.00),
            IntegrationPointType( std::sqrt((7.00-2.00*std::sqrt(7)) / 21.00), (14.00+std::sqrt(7)) / 30.00),
            IntegrationPointType( std::sqrt((7.00+2.00*std::sqrt(7)) / 21.00), (14.00-std::sqrt(7)) / 30.00),
            IntegrationPointType( 1.00, 1.00 / 15.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 6 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints6



class SolidMechanicsLineGaussLobattoIntegrationPoints7
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints7);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 7> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 7;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00, 1.00 / 21.00),
            IntegrationPointType(-std::sqrt((5.00/11.00) + (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 - 7.00*std::sqrt(15.00)) / 350.00),
            IntegrationPointType(-std::sqrt((5.00/11.00) - (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 + 7.00*std::sqrt(15.00)) / 350.00),
            IntegrationPointType(0.00, 256.00/525.00),
            IntegrationPointType(std::sqrt((5.00/11.00) - (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 + 7.00*std::sqrt(15.00)) / 350.00),
            IntegrationPointType(std::sqrt((5.00/11.00) + (2.00/11.00)*std::sqrt(5.00/3.00)), (124.00 - 7.00*std::sqrt(15.00)) / 350.00),
            IntegrationPointType(-1.00, 1.00 / 21.00)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 7 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints7



class SolidMechanicsLineGaussLobattoIntegrationPoints8
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints8);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 8> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 8;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00             , 0.035714285714286),
            IntegrationPointType(-0.871740148509607, 0.210704227143506),
            IntegrationPointType(-0.591700181433142, 0.341122692483504),
            IntegrationPointType(-0.209299217902479, 0.412458794658704),
            IntegrationPointType( 0.209299217902479, 0.412458794658704),
            IntegrationPointType( 0.591700181433142, 0.341122692483504),
            IntegrationPointType( 0.871740148509607, 0.210704227143506),
            IntegrationPointType( 1.00             , 0.035714285714286)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 8 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints8



class SolidMechanicsLineGaussLobattoIntegrationPoints9
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints9);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 9> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 9;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00             , 0.027777777777778),
            IntegrationPointType(-0.899757995411460, 0.165495361560806),
            IntegrationPointType(-0.677186279510738, 0.274538712500162),
            IntegrationPointType(-0.363117463826178, 0.346428510973046),
            IntegrationPointType( 0.00             , 0.371519274376417),
            IntegrationPointType( 0.363117463826178, 0.346428510973046),
            IntegrationPointType( 0.677186279510738, 0.274538712500162),
            IntegrationPointType( 0.899757995411460, 0.165495361560806),
            IntegrationPointType( 1.00             , 0.027777777777778)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 9 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints9



class SolidMechanicsLineGaussLobattoIntegrationPoints10
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SolidMechanicsLineGaussLobattoIntegrationPoints10);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 1;

    typedef IntegrationPoint<1> IntegrationPointType;

    typedef std::array<IntegrationPointType, 10> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 10;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType(-1.00             , 0.022222222222222),
            IntegrationPointType(-0.919533908166459, 0.133305990851070),
            IntegrationPointType(-0.738773865105505, 0.224889342063126),
            IntegrationPointType(-0.477924949810444, 0.292042683679684),
            IntegrationPointType(-0.165278957666387, 0.327539761183897),
            IntegrationPointType( 0.165278957666387, 0.327539761183897),
            IntegrationPointType( 0.477924949810444, 0.292042683679684),
            IntegrationPointType( 0.738773865105505, 0.224889342063126),
            IntegrationPointType( 0.919533908166459, 0.133305990851070),
            IntegrationPointType( 1.00             , 0.022222222222222)
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Line Gauss-Lobatto quadrature 10 ";
        return buffer.str();
    }


}; // Class SolidMechanicsLineGaussLobattoIntegrationPoints10


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.



