//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Stefano Zaghi - Massimo Petracca $
//   Date:                $Date: 2007-03-19 10:49:46 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(INTERFACE_INTEGRATION_POINTS_H_INCLUDED )
#define  INTERFACE_INTEGRATION_POINTS_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <numeric>
#include <cstddef>


// External includes
#include <boost/array.hpp>


// Project includes
#include "includes/define.h"
#include "integration/integration_point.h"


namespace Kratos
{


class KRATOS_API(KRATOS_CORE) QuadrilateralInterfaceLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralInterfaceLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 4> IntegrationPointsArrayType;

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
        buffer << "Quadrilateral interface integration 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;
};



class KRATOS_API(KRATOS_CORE) HexaedralInterfaceLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexaedralInterfaceLobattoIntegrationPoints2);
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
        buffer << "Hexahedral interface quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexaedralInterfaceLobattoIntegrationPoints2



class KRATOS_API(KRATOS_CORE) PrismInterfaceLobattoIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrismInterfaceLobattoIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef boost::array<IntegrationPointType, 6> IntegrationPointsArrayType;

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
        buffer << "Prism interface quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class PrismInterfaceLobattoIntegrationPoints2



}

#endif // INTERFACE_INTEGRATION_POINTS_H_INCLUDED  defined 


