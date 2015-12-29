//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   JMCarbonell $
//   Date:                $Date:   December 2015 $
//   Revision:            $Revision:         1.7 $
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


class KRATOS_API(KRATOS_CORE) HexahedronGaussLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 8> IntegrationPointsArrayType;

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


