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
class KRATOS_API(KRATOS_CORE) HexahedronGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( -1.00 , -1.00, 0.00, 1.00 );
        msIntegrationPoints[1] = IntegrationPointType(  1.00 , -1.00, 0.00, 1.00 );
        msIntegrationPoints[2] = IntegrationPointType(  1.00 ,  1.00, 0.00, 1.00 );
        msIntegrationPoints[3] = IntegrationPointType( -1.00 ,  1.00, 0.00, 1.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLobattoIntegrationPoints1

class KRATOS_API(KRATOS_CORE) HexahedronGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 8> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 8;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        // This is added to solve the problem of static initialization. Pooyan.
        msIntegrationPoints[0] = IntegrationPointType( -1.00 , -1.00, -1.00, 0.50 );
        msIntegrationPoints[1] = IntegrationPointType(  1.00 , -1.00, -1.00, 0.50 );
        msIntegrationPoints[2] = IntegrationPointType(  1.00 ,  1.00, -1.00, 0.50 );
        msIntegrationPoints[3] = IntegrationPointType( -1.00 ,  1.00, -1.00, 0.50 );

        msIntegrationPoints[4] = IntegrationPointType( -1.00 , -1.00,  1.00, 0.50 );
        msIntegrationPoints[5] = IntegrationPointType(  1.00 , -1.00,  1.00, 0.50 );
        msIntegrationPoints[6] = IntegrationPointType(  1.00 ,  1.00,  1.00, 0.50 );
        msIntegrationPoints[7] = IntegrationPointType( -1.00 ,  1.00,  1.00, 0.50 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Lobatto quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLobattoIntegrationPoints2



}

#endif // KRATOS_HEXAHEDRON_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED defined


