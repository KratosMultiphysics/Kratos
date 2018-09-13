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

#if !defined(KRATOS_QUADRILATERAL_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

//TO BE COMPLETED: Only the needed ones have been implemented

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) QuadrilateralGaussLobattoIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 2> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 2;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( -1.00 , 0.00, 1.00 );
        msIntegrationPoints[1] = IntegrationPointType(  1.00 , 0.00, 1.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class QuadrilateralGaussLobattoIntegrationPoints1

class KRATOS_API(KRATOS_CORE) QuadrilateralGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralGaussLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 4> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 4;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        // This is added to solve the problem of static initialization. Pooyan.
        msIntegrationPoints[0] = IntegrationPointType( -1.00 , -1.00, 0.50 );
        msIntegrationPoints[1] = IntegrationPointType(  1.00 , -1.00, 0.50 );
        msIntegrationPoints[2] = IntegrationPointType(  1.00 ,  1.00, 0.50 );
        msIntegrationPoints[3] = IntegrationPointType( -1.00 ,  1.00, 0.50 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Gauss-Lobatto integration 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;
}; // Class QuadrilateralGaussLobattoIntegrationPoints2



}

#endif // KRATOS_QUADRILATERAL_GAUSS_LOBATTO_INTEGRATION_POINTS_H_INCLUDED defined


