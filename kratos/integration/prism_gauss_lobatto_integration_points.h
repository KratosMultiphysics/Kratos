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

#if !defined(KRATOS_PRISM_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_PRISM_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

//TO BE COMPLETED: Only the needed ones have been implemented

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) PrismGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.0 , 0.0 , 0.5 , 1.00 / 6.00 );
        msIntegrationPoints[1] = IntegrationPointType( 1.0 , 0.0 , 0.5 , 1.00 / 6.00 );
        msIntegrationPoints[2] = IntegrationPointType( 0.0 , 1.0 , 0.5 , 1.00 / 6.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Lobatto quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLobattoIntegrationPoints1

class KRATOS_API(KRATOS_CORE) PrismGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismGaussLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        // This is added to solve the problem of static initialization. Pooyan.
        msIntegrationPoints[0] = IntegrationPointType( 0.0 , 0.0 , 0.0 , 1.00 / 12.00 );
        msIntegrationPoints[1] = IntegrationPointType( 1.0 , 0.0 , 0.0 , 1.00 / 12.00 );
        msIntegrationPoints[2] = IntegrationPointType( 0.0 , 1.0 , 0.0 , 1.00 / 12.00 );

        msIntegrationPoints[3] = IntegrationPointType( 0.0 , 0.0 , 1.0 , 1.00 / 12.00 );
        msIntegrationPoints[4] = IntegrationPointType( 1.0 , 0.0 , 1.0 , 1.00 / 12.00 );
        msIntegrationPoints[5] = IntegrationPointType( 0.0 , 1.0 , 1.0 , 1.00 / 12.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Prism Gauss-Lobatto quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismGaussLobattoIntegrationPoints2



}

#endif // KRATOS_PRISM_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED defined


