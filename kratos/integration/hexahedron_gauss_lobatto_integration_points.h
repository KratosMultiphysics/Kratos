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

#if !defined(KRATOS_HEXAHEDRON_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_HEXAHEDRON_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

//TO BE COMPLETED: Only the needed ones have been implemented

namespace Kratos
{

//TODO
class HexahedronGaussLobattoIntegrationPoints0
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLobattoIntegrationPoints0);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( -1.0 , -1.0, 0.0, 1.0 ),
            IntegrationPointType(  1.0 , -1.0, 0.0, 1.0 ),
            IntegrationPointType(  1.0 ,  1.0, 0.0, 1.0 ),
            IntegrationPointType( -1.0 ,  1.0, 0.0, 1.0 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Lobatto quadrature 0 ";
        return buffer.str();
    }
}; // Class HexahedronGaussLobattoIntegrationPoints0

class HexahedronGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 8> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 8;
    }

    static const IntegrationPointsArrayType& IntegrationPoints()
    {
        static const IntegrationPointsArrayType s_integration_points{{
            IntegrationPointType( -1.0 , -1.0, -1.0, 1.0 ),
            IntegrationPointType(  1.0 , -1.0, -1.0, 1.0 ),
            IntegrationPointType(  1.0 ,  1.0, -1.0, 1.0 ),
            IntegrationPointType( -1.0 ,  1.0, -1.0, 1.0 ),
            IntegrationPointType( -1.0 , -1.0,  1.0, 1.0 ),
            IntegrationPointType(  1.0 , -1.0,  1.0, 1.0 ),
            IntegrationPointType(  1.0 ,  1.0,  1.0, 1.0 ),
            IntegrationPointType( -1.0 ,  1.0,  1.0, 1.0 )
        }};
        return s_integration_points;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }
}; // Class HexahedronGaussLobattoIntegrationPoints1



}

#endif // KRATOS_HEXAHEDRON_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED defined


