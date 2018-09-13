//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_QUADRILATERAL_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_QUADRILATERAL_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{

class KRATOS_API(KRATOS_CORE) QuadrilateralCollocationIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints1);
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
        msIntegrationPoints[0]=IntegrationPointType(-0.500000000000,-0.500000000000,1.00000000000);
        msIntegrationPoints[1]=IntegrationPointType(-0.500000000000,0.500000000000,1.00000000000);
        msIntegrationPoints[2]=IntegrationPointType(0.500000000000,-0.500000000000,1.00000000000);
        msIntegrationPoints[3]=IntegrationPointType(0.500000000000,0.500000000000,1.00000000000);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class QuadrilateralCollocationIntegrationPoints1

class KRATOS_API(KRATOS_CORE) QuadrilateralCollocationIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints2);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 9> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 9;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(-0.666666666667,-0.666666666667,0.444444444444);
        msIntegrationPoints[1]=IntegrationPointType(-0.666666666667,0,0.444444444444);
        msIntegrationPoints[2]=IntegrationPointType(-0.666666666667,0.666666666667,0.444444444444);
        msIntegrationPoints[3]=IntegrationPointType(0,-0.666666666667,0.444444444444);
        msIntegrationPoints[4]=IntegrationPointType(0,0,0.444444444444);
        msIntegrationPoints[5]=IntegrationPointType(0,0.666666666667,0.444444444444);
        msIntegrationPoints[6]=IntegrationPointType(0.666666666667,-0.666666666667,0.444444444444);
        msIntegrationPoints[7]=IntegrationPointType(0.666666666667,0,0.444444444444);
        msIntegrationPoints[8]=IntegrationPointType(0.666666666667,0.666666666667,0.444444444444);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class QuadrilateralCollocationIntegrationPoints2

class KRATOS_API(KRATOS_CORE) QuadrilateralCollocationIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 16> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 16;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(-0.750000000000,-0.750000000000,0.250000000000);
        msIntegrationPoints[1]=IntegrationPointType(-0.750000000000,-0.250000000000,0.250000000000);
        msIntegrationPoints[2]=IntegrationPointType(-0.750000000000,0.250000000000,0.250000000000);
        msIntegrationPoints[3]=IntegrationPointType(-0.750000000000,0.750000000000,0.250000000000);
        msIntegrationPoints[4]=IntegrationPointType(-0.250000000000,-0.750000000000,0.250000000000);
        msIntegrationPoints[5]=IntegrationPointType(-0.250000000000,-0.250000000000,0.250000000000);
        msIntegrationPoints[6]=IntegrationPointType(-0.250000000000,0.250000000000,0.250000000000);
        msIntegrationPoints[7]=IntegrationPointType(-0.250000000000,0.750000000000,0.250000000000);
        msIntegrationPoints[8]=IntegrationPointType(0.250000000000,-0.750000000000,0.250000000000);
        msIntegrationPoints[9]=IntegrationPointType(0.250000000000,-0.250000000000,0.250000000000);
        msIntegrationPoints[10]=IntegrationPointType(0.250000000000,0.250000000000,0.250000000000);
        msIntegrationPoints[11]=IntegrationPointType(0.250000000000,0.750000000000,0.250000000000);
        msIntegrationPoints[12]=IntegrationPointType(0.750000000000,-0.750000000000,0.250000000000);
        msIntegrationPoints[13]=IntegrationPointType(0.750000000000,-0.250000000000,0.250000000000);
        msIntegrationPoints[14]=IntegrationPointType(0.750000000000,0.250000000000,0.250000000000);
        msIntegrationPoints[15]=IntegrationPointType(0.750000000000,0.750000000000,0.250000000000);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class QuadrilateralCollocationIntegrationPoints3

class KRATOS_API(KRATOS_CORE) QuadrilateralCollocationIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 25> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 25;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0]=IntegrationPointType(-0.800000000000,-0.800000000000,0.160000000000);
        msIntegrationPoints[1]=IntegrationPointType(-0.800000000000,-0.400000000000,0.160000000000);
        msIntegrationPoints[2]=IntegrationPointType(-0.800000000000, 0.000000000000,0.160000000000);
        msIntegrationPoints[3]=IntegrationPointType(-0.800000000000,0.400000000000,0.160000000000);
        msIntegrationPoints[4]=IntegrationPointType(-0.800000000000,0.800000000000,0.160000000000);
        msIntegrationPoints[5]=IntegrationPointType(-0.400000000000,-0.800000000000,0.160000000000);
        msIntegrationPoints[6]=IntegrationPointType(-0.400000000000,-0.400000000000,0.160000000000);
        msIntegrationPoints[7]=IntegrationPointType(-0.400000000000,0.000000000000,0.160000000000);
        msIntegrationPoints[8]=IntegrationPointType(-0.400000000000,0.400000000000,0.160000000000);
        msIntegrationPoints[9]=IntegrationPointType(-0.400000000000,0.800000000000,0.160000000000);
        msIntegrationPoints[10]=IntegrationPointType(0.0000000000000,-0.800000000000,0.160000000000);
        msIntegrationPoints[11]=IntegrationPointType(0.0000000000000,-0.400000000000,0.160000000000);
        msIntegrationPoints[12]=IntegrationPointType(0.0000000000000, 0.000000000000,0.160000000000);
        msIntegrationPoints[13]=IntegrationPointType(0.0000000000000, 0.400000000000,0.160000000000);
        msIntegrationPoints[14]=IntegrationPointType(0.0000000000000,0.800000000000,0.160000000000);
        msIntegrationPoints[15]=IntegrationPointType(0.400000000000,-0.800000000000,0.160000000000);
        msIntegrationPoints[16]=IntegrationPointType(0.400000000000,-0.400000000000,0.160000000000);
        msIntegrationPoints[17]=IntegrationPointType(0.400000000000, 0.000000000000,0.160000000000);
        msIntegrationPoints[18]=IntegrationPointType(0.400000000000,0.400000000000,0.160000000000);
        msIntegrationPoints[19]=IntegrationPointType(0.400000000000,0.800000000000,0.160000000000);
        msIntegrationPoints[20]=IntegrationPointType(0.800000000000,-0.800000000000,0.160000000000);
        msIntegrationPoints[21]=IntegrationPointType(0.800000000000,-0.400000000000,0.160000000000);
        msIntegrationPoints[22]=IntegrationPointType(0.800000000000, 0.000000000000,0.160000000000);
        msIntegrationPoints[23]=IntegrationPointType(0.800000000000,0.400000000000,0.160000000000);
        msIntegrationPoints[24]=IntegrationPointType(0.800000000000,0.800000000000,0.160000000000);
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Quadrilateral Collocation quadrature 4 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class QuadrilateralCollocationIntegrationPoints4

class KRATOS_API(KRATOS_CORE) QuadrilateralCollocationIntegrationPoints5 {
public:
	KRATOS_CLASS_POINTER_DEFINITION(QuadrilateralCollocationIntegrationPoints5);
	typedef std::size_t SizeType;

	static const unsigned int Dimension = 2;

	typedef IntegrationPoint<2> IntegrationPointType;

	typedef std::array<IntegrationPointType, 36> IntegrationPointsArrayType;

	typedef IntegrationPointType::PointType PointType;

	static SizeType IntegrationPointsNumber()
        {
            return 36;
        }

	static IntegrationPointsArrayType& IntegrationPoints()
	{
            msIntegrationPoints[0]=IntegrationPointType(-0.833333333333,-0.833333333333,0.111111111111);
            msIntegrationPoints[1]=IntegrationPointType(-0.833333333333,-0.500000000000,0.111111111111);
            msIntegrationPoints[2]=IntegrationPointType(-0.833333333333,-0.166666666667,0.111111111111);
            msIntegrationPoints[3]=IntegrationPointType(-0.833333333333,0.166666666667,0.111111111111);
            msIntegrationPoints[4]=IntegrationPointType(-0.833333333333,0.500000000000,0.111111111111);
            msIntegrationPoints[5]=IntegrationPointType(-0.833333333333,0.833333333333,0.111111111111);
            msIntegrationPoints[6]=IntegrationPointType(-0.500000000000,-0.833333333333,0.111111111111);
            msIntegrationPoints[7]=IntegrationPointType(-0.500000000000,-0.500000000000,0.111111111111);
            msIntegrationPoints[8]=IntegrationPointType(-0.500000000000,-0.166666666667,0.111111111111);
            msIntegrationPoints[9]=IntegrationPointType(-0.500000000000,0.166666666667,0.111111111111);
            msIntegrationPoints[10]=IntegrationPointType(-0.500000000000,0.500000000000,0.111111111111);
            msIntegrationPoints[11]=IntegrationPointType(-0.500000000000,0.833333333333,0.111111111111);
            msIntegrationPoints[12]=IntegrationPointType(-0.166666666667,-0.833333333333,0.111111111111);
            msIntegrationPoints[13]=IntegrationPointType(-0.166666666667,-0.500000000000,0.111111111111);
            msIntegrationPoints[14]=IntegrationPointType(-0.166666666667,-0.166666666667,0.111111111111);
            msIntegrationPoints[15]=IntegrationPointType(-0.166666666667,0.166666666667,0.111111111111);
            msIntegrationPoints[16]=IntegrationPointType(-0.166666666667,0.500000000000,0.111111111111);
            msIntegrationPoints[17]=IntegrationPointType(-0.166666666667,0.833333333333,0.111111111111);
            msIntegrationPoints[18]=IntegrationPointType(0.166666666667,-0.833333333333,0.111111111111);
            msIntegrationPoints[19]=IntegrationPointType(0.166666666667,-0.500000000000,0.111111111111);
            msIntegrationPoints[20]=IntegrationPointType(0.166666666667,-0.166666666667,0.111111111111);
            msIntegrationPoints[21]=IntegrationPointType(0.166666666667,0.166666666667,0.111111111111);
            msIntegrationPoints[22]=IntegrationPointType(0.166666666667,0.500000000000,0.111111111111);
            msIntegrationPoints[23]=IntegrationPointType(0.166666666667,0.833333333333,0.111111111111);
            msIntegrationPoints[24]=IntegrationPointType(0.500000000000,-0.833333333333,0.111111111111);
            msIntegrationPoints[25]=IntegrationPointType(0.500000000000,-0.500000000000,0.111111111111);
            msIntegrationPoints[26]=IntegrationPointType(0.500000000000,-0.166666666667,0.111111111111);
            msIntegrationPoints[27]=IntegrationPointType(0.500000000000,0.166666666667,0.111111111111);
            msIntegrationPoints[28]=IntegrationPointType(0.500000000000,0.500000000000,0.111111111111);
            msIntegrationPoints[29]=IntegrationPointType(0.500000000000,0.833333333333,0.111111111111);
            msIntegrationPoints[30]=IntegrationPointType(0.833333333333,-0.833333333333,0.111111111111);
            msIntegrationPoints[31]=IntegrationPointType(0.833333333333,-0.500000000000,0.111111111111);
            msIntegrationPoints[32]=IntegrationPointType(0.833333333333,-0.166666666667,0.111111111111);
            msIntegrationPoints[33]=IntegrationPointType(0.833333333333,0.166666666667,0.111111111111);
            msIntegrationPoints[34]=IntegrationPointType(0.833333333333,0.500000000000,0.111111111111);
            msIntegrationPoints[35]=IntegrationPointType(0.833333333333,0.833333333333,0.111111111111);
            return msIntegrationPoints;
	}

	std::string Info() const
	{
		std::stringstream buffer;
		buffer << "Quadrilateral Collocation quadrature 5 ";
		return buffer.str();
	}
protected:

private:

	static IntegrationPointsArrayType msIntegrationPoints;

}; // Class QuadrilateralCollocationIntegrationPoints5


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_QUADRILATERAL_COLLOCATION_INTEGRATION_POINTS_H_INCLUDED  defined


