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


#if !defined(KRATOS_HEXAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_HEXAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 1> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
       msIntegrationPoints[0] = IntegrationPointType( 0.00 , 0.00, 0.00 , 8.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints1

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints2);
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
       msIntegrationPoints[0] = IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
       msIntegrationPoints[1] = IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
       msIntegrationPoints[2] = IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
       msIntegrationPoints[3] = IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0), -1.00/std::sqrt(3.0), 1.00 );
       msIntegrationPoints[4] = IntegrationPointType( -1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
       msIntegrationPoints[5] = IntegrationPointType(  1.00/std::sqrt(3.0) , -1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
       msIntegrationPoints[6] = IntegrationPointType(  1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );
       msIntegrationPoints[7] = IntegrationPointType( -1.00/std::sqrt(3.0) ,  1.00/std::sqrt(3.0),  1.00/std::sqrt(3.0), 1.00 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexahedron Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints2

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 27> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 27;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
       msIntegrationPoints[ 0] = IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );
       msIntegrationPoints[ 1] = IntegrationPointType(                   0.0 , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 200.00/729.00 );
       msIntegrationPoints[ 2] = IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );

       msIntegrationPoints[ 3] = IntegrationPointType( -std::sqrt(3.00/5.00) ,                   0.0, -std::sqrt(3.00/5.00), 200.00/729.00 );
       msIntegrationPoints[ 4] = IntegrationPointType(                   0.0 ,                   0.0, -std::sqrt(3.00/5.00), 320.00/729.00 );
       msIntegrationPoints[ 5] = IntegrationPointType(  std::sqrt(3.00/5.00) ,                   0.0, -std::sqrt(3.00/5.00), 200.00/729.00 );

       msIntegrationPoints[ 6] = IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );
       msIntegrationPoints[ 7] = IntegrationPointType(                   0.0 ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 200.00/729.00 );
       msIntegrationPoints[ 8] = IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00), -std::sqrt(3.00/5.00), 125.00/729.00 );

       msIntegrationPoints[ 9] = IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );
       msIntegrationPoints[10] = IntegrationPointType(                   0.0 , -std::sqrt(3.00/5.00),                   0.0, 320.00/729.00 );
       msIntegrationPoints[11] = IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );

       msIntegrationPoints[12] = IntegrationPointType( -std::sqrt(3.00/5.00) ,                   0.0,                   0.0, 320.00/729.00 );
       msIntegrationPoints[13] = IntegrationPointType(                   0.0 ,                   0.0,                   0.0, 512.00/729.00 );
       msIntegrationPoints[14] = IntegrationPointType(  std::sqrt(3.00/5.00) ,                   0.0,                   0.0, 320.00/729.00 );

       msIntegrationPoints[15] = IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );
       msIntegrationPoints[16] = IntegrationPointType(                   0.0 ,  std::sqrt(3.00/5.00),                   0.0, 320.00/729.00 );
       msIntegrationPoints[17] = IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),                   0.0, 200.00/729.00 );

       msIntegrationPoints[18] = IntegrationPointType( -std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );
       msIntegrationPoints[19] = IntegrationPointType(                   0.0 , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 200.00/729.00 );
       msIntegrationPoints[20] = IntegrationPointType(  std::sqrt(3.00/5.00) , -std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );

       msIntegrationPoints[21] = IntegrationPointType( -std::sqrt(3.00/5.00) ,                   0.0,  std::sqrt(3.00/5.00), 200.00/729.00 );
       msIntegrationPoints[22] = IntegrationPointType(                   0.0 ,                   0.0,  std::sqrt(3.00/5.00), 320.00/729.00 );
       msIntegrationPoints[23] = IntegrationPointType(  std::sqrt(3.00/5.00) ,                   0.0,  std::sqrt(3.00/5.00), 200.00/729.00 );

       msIntegrationPoints[24] = IntegrationPointType( -std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );
       msIntegrationPoints[25] = IntegrationPointType(                   0.0 ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 200.00/729.00 );
       msIntegrationPoints[26] = IntegrationPointType(  std::sqrt(3.00/5.00) ,  std::sqrt(3.00/5.00),  std::sqrt(3.00/5.00), 125.00/729.00 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Hexadra Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints3

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints4
{
public:
   KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints4);
   typedef std::size_t SizeType;

   static const unsigned int Dimension = 3;

   typedef IntegrationPoint<3> IntegrationPointType;

   typedef std::array<IntegrationPointType, 64> IntegrationPointsArrayType;

   typedef IntegrationPointType::PointType PointType;

   static SizeType IntegrationPointsNumber()
   {
       return 64;
   }

   static IntegrationPointsArrayType& IntegrationPoints()
   {
       msIntegrationPoints[0] = IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454);
       msIntegrationPoints[1] = IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[2] = IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[3] = IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454);
       msIntegrationPoints[4] = IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[5] = IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[6] = IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[7] = IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[8] = IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[9] = IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[10] = IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , -0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[11] = IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[12] = IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454);
       msIntegrationPoints[13] = IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[14] = IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , -0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[15] = IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , -0.86113631159405257522 , 0.04209147749053145454);
       msIntegrationPoints[16] = IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[17] = IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[18] = IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[19] = IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[20] = IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[21] = IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[22] = IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[23] = IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[24] = IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[25] = IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[26] = IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , -0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[27] = IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[28] = IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[29] = IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[30] = IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , -0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[31] = IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , -0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[32] = IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[33] = IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[34] = IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[35] = IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[36] = IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[37] = IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[38] = IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[39] = IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[40] = IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[41] = IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[42] = IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , 0.33998104358485626480 , 0.27735296695391298990);
       msIntegrationPoints[43] = IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[44] = IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[45] = IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[46] = IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , 0.33998104358485626480 , 0.14794033605678130087);
       msIntegrationPoints[47] = IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , 0.33998104358485626480 , 0.07891151579507055098);
       msIntegrationPoints[48] = IntegrationPointType(-0.86113631159405257522 , -0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454);
       msIntegrationPoints[49] = IntegrationPointType(-0.33998104358485626480 , -0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[50] = IntegrationPointType(0.33998104358485626480 , -0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[51] = IntegrationPointType(0.86113631159405257522 , -0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454);
       msIntegrationPoints[52] = IntegrationPointType(-0.86113631159405257522 , -0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[53] = IntegrationPointType(-0.33998104358485626480 , -0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[54] = IntegrationPointType(0.33998104358485626480 , -0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[55] = IntegrationPointType(0.86113631159405257522 , -0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[56] = IntegrationPointType(-0.86113631159405257522 , 0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[57] = IntegrationPointType(-0.33998104358485626480 , 0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[58] = IntegrationPointType(0.33998104358485626480 , 0.33998104358485626480 , 0.86113631159405257522 , 0.14794033605678130087);
       msIntegrationPoints[59] = IntegrationPointType(0.86113631159405257522 , 0.33998104358485626480 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[60] = IntegrationPointType(-0.86113631159405257522 , 0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454);
       msIntegrationPoints[61] = IntegrationPointType(-0.33998104358485626480 , 0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[62] = IntegrationPointType(0.33998104358485626480 , 0.86113631159405257522 , 0.86113631159405257522 , 0.07891151579507055098);
       msIntegrationPoints[63] = IntegrationPointType(0.86113631159405257522 , 0.86113631159405257522 , 0.86113631159405257522 , 0.04209147749053145454);

       return msIntegrationPoints;
   }

   std::string Info() const
   {
       std::stringstream buffer;
       buffer << "Hexadra Gauss-Legendre quadrature 4 ";
       return buffer.str();
   }
protected:

private:

   static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints4

class KRATOS_API(KRATOS_CORE) HexahedronGaussLegendreIntegrationPoints5
{
public:
   KRATOS_CLASS_POINTER_DEFINITION(HexahedronGaussLegendreIntegrationPoints5);
   typedef std::size_t SizeType;

   static const unsigned int Dimension = 3;

   typedef IntegrationPoint<3> IntegrationPointType;

   typedef std::array<IntegrationPointType, 125> IntegrationPointsArrayType;

   typedef IntegrationPointType::PointType PointType;

   static SizeType IntegrationPointsNumber()
   {
       return 125;
   }

   static IntegrationPointsArrayType& IntegrationPoints()
   {
       msIntegrationPoints[0] = IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092);
       msIntegrationPoints[1] = IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[2] = IntegrationPointType(0 , -0.90617984593866399280 , -0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[3] = IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[4] = IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092);
       msIntegrationPoints[5] = IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[6] = IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[7] = IntegrationPointType(0 , -0.53846931010568309104 , -0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[8] = IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[9] = IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[10] = IntegrationPointType(-0.90617984593866399280 , 0 , -0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[11] = IntegrationPointType(-0.53846931010568309104 , 0 , -0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[12] = IntegrationPointType(0 , 0 , -0.90617984593866399280 , 0.07667773006934522489);
       msIntegrationPoints[13] = IntegrationPointType(0.53846931010568309104 , 0 , -0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[14] = IntegrationPointType(0.90617984593866399280 , 0 , -0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[15] = IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[16] = IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[17] = IntegrationPointType(0 , 0.53846931010568309104 , -0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[18] = IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , -0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[19] = IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[20] = IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092);
       msIntegrationPoints[21] = IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[22] = IntegrationPointType(0 , 0.90617984593866399280 , -0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[23] = IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , -0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[24] = IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , -0.90617984593866399280 , 0.013299736420632648092);
       msIntegrationPoints[25] = IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[26] = IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[27] = IntegrationPointType(0 , -0.90617984593866399280 , -0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[28] = IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[29] = IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[30] = IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[31] = IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[32] = IntegrationPointType(0 , -0.53846931010568309104 , -0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[33] = IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[34] = IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[35] = IntegrationPointType(-0.90617984593866399280 , 0 , -0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[36] = IntegrationPointType(-0.53846931010568309104 , 0 , -0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[37] = IntegrationPointType(0 , 0 , -0.53846931010568309104 , 0.15490078296220484370);
       msIntegrationPoints[38] = IntegrationPointType(0.53846931010568309104 , 0 , -0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[39] = IntegrationPointType(0.90617984593866399280 , 0 , -0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[40] = IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[41] = IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[42] = IntegrationPointType(0 , 0.53846931010568309104 , -0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[43] = IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , -0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[44] = IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[45] = IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[46] = IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[47] = IntegrationPointType(0 , 0.90617984593866399280 , -0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[48] = IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , -0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[49] = IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , -0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[50] = IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , 0 , 0.031934207352848290676);
       msIntegrationPoints[51] = IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , 0 , 0.06451200000000000000);
       msIntegrationPoints[52] = IntegrationPointType(0 , -0.90617984593866399280 , 0 , 0.07667773006934522489);
       msIntegrationPoints[53] = IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , 0 , 0.06451200000000000000);
       msIntegrationPoints[54] = IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , 0 , 0.031934207352848290676);
       msIntegrationPoints[55] = IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , 0 , 0.06451200000000000000);
       msIntegrationPoints[56] = IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , 0 , 0.13032414106964827997);
       msIntegrationPoints[57] = IntegrationPointType(0 , -0.53846931010568309104 , 0 , 0.15490078296220484370);
       msIntegrationPoints[58] = IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , 0 , 0.13032414106964827997);
       msIntegrationPoints[59] = IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , 0 , 0.06451200000000000000);
       msIntegrationPoints[60] = IntegrationPointType(-0.90617984593866399280 , 0 , 0 , 0.07667773006934522489);
       msIntegrationPoints[61] = IntegrationPointType(-0.53846931010568309104 , 0 , 0 , 0.15490078296220484370);
       msIntegrationPoints[62] = IntegrationPointType(0 , 0 , 0 , 0.18411210973936899863);
       msIntegrationPoints[63] = IntegrationPointType(0.53846931010568309104 , 0 , 0 , 0.15490078296220484370);
       msIntegrationPoints[64] = IntegrationPointType(0.90617984593866399280 , 0 , 0 , 0.07667773006934522489);
       msIntegrationPoints[65] = IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , 0 , 0.06451200000000000000);
       msIntegrationPoints[66] = IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , 0 , 0.13032414106964827997);
       msIntegrationPoints[67] = IntegrationPointType(0 , 0.53846931010568309104 , 0 , 0.15490078296220484370);
       msIntegrationPoints[68] = IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , 0 , 0.13032414106964827997);
       msIntegrationPoints[69] = IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , 0 , 0.06451200000000000000);
       msIntegrationPoints[70] = IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , 0 , 0.031934207352848290676);
       msIntegrationPoints[71] = IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , 0 , 0.06451200000000000000);
       msIntegrationPoints[72] = IntegrationPointType(0 , 0.90617984593866399280 , 0 , 0.07667773006934522489);
       msIntegrationPoints[73] = IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , 0 , 0.06451200000000000000);
       msIntegrationPoints[74] = IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , 0 , 0.031934207352848290676);
       msIntegrationPoints[75] = IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[76] = IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[77] = IntegrationPointType(0 , -0.90617984593866399280 , 0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[78] = IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[79] = IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[80] = IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[81] = IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[82] = IntegrationPointType(0 , -0.53846931010568309104 , 0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[83] = IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[84] = IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[85] = IntegrationPointType(-0.90617984593866399280 , 0 , 0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[86] = IntegrationPointType(-0.53846931010568309104 , 0 , 0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[87] = IntegrationPointType(0 , 0 , 0.53846931010568309104 , 0.15490078296220484370);
       msIntegrationPoints[88] = IntegrationPointType(0.53846931010568309104 , 0 , 0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[89] = IntegrationPointType(0.90617984593866399280 , 0 , 0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[90] = IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[91] = IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[92] = IntegrationPointType(0 , 0.53846931010568309104 , 0.53846931010568309104 , 0.13032414106964827997);
       msIntegrationPoints[93] = IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , 0.53846931010568309104 , 0.10964684245453881967);
       msIntegrationPoints[94] = IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[95] = IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[96] = IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[97] = IntegrationPointType(0 , 0.90617984593866399280 , 0.53846931010568309104 , 0.06451200000000000000);
       msIntegrationPoints[98] = IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , 0.53846931010568309104 , 0.05427649123462815748);
       msIntegrationPoints[99] = IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , 0.53846931010568309104 , 0.026867508765371842524);
       msIntegrationPoints[100] = IntegrationPointType(-0.90617984593866399280 , -0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092);
       msIntegrationPoints[101] = IntegrationPointType(-0.53846931010568309104 , -0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[102] = IntegrationPointType(0 , -0.90617984593866399280 , 0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[103] = IntegrationPointType(0.53846931010568309104 , -0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[104] = IntegrationPointType(0.90617984593866399280 , -0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092);
       msIntegrationPoints[105] = IntegrationPointType(-0.90617984593866399280 , -0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[106] = IntegrationPointType(-0.53846931010568309104 , -0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[107] = IntegrationPointType(0 , -0.53846931010568309104 , 0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[108] = IntegrationPointType(0.53846931010568309104 , -0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[109] = IntegrationPointType(0.90617984593866399280 , -0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[110] = IntegrationPointType(-0.90617984593866399280 , 0 , 0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[111] = IntegrationPointType(-0.53846931010568309104 , 0 , 0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[112] = IntegrationPointType(0 , 0 , 0.90617984593866399280 , 0.07667773006934522489);
       msIntegrationPoints[113] = IntegrationPointType(0.53846931010568309104 , 0 , 0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[114] = IntegrationPointType(0.90617984593866399280 , 0 , 0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[115] = IntegrationPointType(-0.90617984593866399280 , 0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[116] = IntegrationPointType(-0.53846931010568309104 , 0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[117] = IntegrationPointType(0 , 0.53846931010568309104 , 0.90617984593866399280 , 0.06451200000000000000);
       msIntegrationPoints[118] = IntegrationPointType(0.53846931010568309104 , 0.53846931010568309104 , 0.90617984593866399280 , 0.05427649123462815748);
       msIntegrationPoints[119] = IntegrationPointType(0.90617984593866399280 , 0.53846931010568309104 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[120] = IntegrationPointType(-0.90617984593866399280 , 0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092);
       msIntegrationPoints[121] = IntegrationPointType(-0.53846931010568309104 , 0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[122] = IntegrationPointType(0 , 0.90617984593866399280 , 0.90617984593866399280 , 0.031934207352848290676);
       msIntegrationPoints[123] = IntegrationPointType(0.53846931010568309104 , 0.90617984593866399280 , 0.90617984593866399280 , 0.026867508765371842524);
       msIntegrationPoints[124] = IntegrationPointType(0.90617984593866399280 , 0.90617984593866399280 , 0.90617984593866399280 , 0.013299736420632648092);

       return msIntegrationPoints;
   }

   std::string Info() const
   {
       std::stringstream buffer;
       buffer << "Hexadra Gauss-Legendre quadrature 5 ";
       return buffer.str();
   }
protected:

private:

   static IntegrationPointsArrayType msIntegrationPoints;

}; // Class HexahedronGaussLegendreIntegrationPoints5

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_HEXAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined


