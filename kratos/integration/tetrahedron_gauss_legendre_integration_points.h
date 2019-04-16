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


#if !defined(KRATOS_TETRAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_TETRAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "integration/quadrature.h"


namespace Kratos
{

class KRATOS_API(KRATOS_CORE) TetrahedronGaussLegendreIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedronGaussLegendreIntegrationPoints1);
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
        msIntegrationPoints[0] = IntegrationPointType( 0.25,0.25,0.25 , 1.00 / 6.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Tetrahedron Gauss-Legendre quadrature 1 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TetrahedronGaussLegendreIntegrationPoints1

class KRATOS_API(KRATOS_CORE) TetrahedronGaussLegendreIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedronGaussLegendreIntegrationPoints2);
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
        msIntegrationPoints[0] = IntegrationPointType( 0.58541020,0.13819660,0.13819660 , 1.00 / 24.00 );
        msIntegrationPoints[1] = IntegrationPointType( 0.13819660,0.58541020,0.13819660 , 1.00 / 24.00 );
        msIntegrationPoints[2] = IntegrationPointType( 0.13819660,0.13819660,0.58541020 , 1.00 / 24.00 );
        msIntegrationPoints[3] = IntegrationPointType( 0.13819660,0.13819660,0.13819660 , 1.00 / 24.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Tetrahedron Gauss-Legendre quadrature 2 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TetrahedronGaussLegendreIntegrationPoints2

class KRATOS_API(KRATOS_CORE) TetrahedronGaussLegendreIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedronGaussLegendreIntegrationPoints3);
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
        msIntegrationPoints[0] = IntegrationPointType( 0.015835909865720057993,0.32805469671142664734,0.32805469671142664734 , 0.02308799441864369039 );
        msIntegrationPoints[1] = IntegrationPointType( 0.32805469671142664734,0.015835909865720057993,0.32805469671142664734 , 0.02308799441864369039 );
        msIntegrationPoints[2] = IntegrationPointType( 0.32805469671142664734,0.32805469671142664734,0.015835909865720057993 , 0.02308799441864369039 );
        msIntegrationPoints[3] = IntegrationPointType( 0.32805469671142664734,0.32805469671142664734,0.32805469671142664734 , 0.02308799441864369039 );
        msIntegrationPoints[4] = IntegrationPointType( 0.67914317820120795168,0.10695227393293068277,0.10695227393293068277 , 0.01857867224802297628 );
        msIntegrationPoints[5] = IntegrationPointType( 0.10695227393293068277,0.67914317820120795168,0.10695227393293068277 , 0.01857867224802297628 );
        msIntegrationPoints[6] = IntegrationPointType( 0.10695227393293068277,0.10695227393293068277,0.67914317820120795168 , 0.01857867224802297628 );
        msIntegrationPoints[7] = IntegrationPointType( 0.10695227393293068277,0.10695227393293068277,0.10695227393293068277 , 0.01857867224802297628 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Tetrahedron Gauss-Legendre quadrature 3 ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TetrahedronGaussLegendreIntegrationPoints3

class KRATOS_API(KRATOS_CORE) TetrahedronGaussLegendreIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedronGaussLegendreIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 14> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 14;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.72179424906732632079,0.092735250310891226402,0.092735250310891226402 , 0.01224884051939365826 );
        msIntegrationPoints[1] = IntegrationPointType( 0.092735250310891226402,0.72179424906732632079,0.092735250310891226402 , 0.01224884051939365826 );
        msIntegrationPoints[2] = IntegrationPointType( 0.092735250310891226402,0.092735250310891226402,0.72179424906732632079 , 0.01224884051939365826 );
        msIntegrationPoints[3] = IntegrationPointType( 0.092735250310891226402,0.092735250310891226402,0.092735250310891226402 , 0.01224884051939365826 );
        msIntegrationPoints[4] = IntegrationPointType( 0.067342242210098170608,0.31088591926330060980,0.31088591926330060980 , 0.01878132095300264180 );
        msIntegrationPoints[5] = IntegrationPointType( 0.31088591926330060980,0.067342242210098170608,0.31088591926330060980 , 0.01878132095300264180 );
        msIntegrationPoints[6] = IntegrationPointType( 0.31088591926330060980,0.31088591926330060980,0.067342242210098170608 , 0.01878132095300264180 );
        msIntegrationPoints[7] = IntegrationPointType( 0.31088591926330060980,0.31088591926330060980,0.31088591926330060980 , 0.01878132095300264180 );
        msIntegrationPoints[8] = IntegrationPointType( 0.045503704125649649492,0.045503704125649649492,0.45449629587435035051 , 0.007091003462846911073 );
        msIntegrationPoints[9] = IntegrationPointType( 0.045503704125649649492,0.45449629587435035051,0.045503704125649649492 , 0.007091003462846911073 );
        msIntegrationPoints[10] = IntegrationPointType( 0.045503704125649649492,0.45449629587435035051,0.45449629587435035051 , 0.007091003462846911073 );
        msIntegrationPoints[11] = IntegrationPointType( 0.45449629587435035051,0.045503704125649649492,0.045503704125649649492 , 0.007091003462846911073 );
        msIntegrationPoints[12] = IntegrationPointType( 0.45449629587435035051,0.045503704125649649492,0.45449629587435035051 , 0.007091003462846911073 );
        msIntegrationPoints[13] = IntegrationPointType( 0.45449629587435035051,0.45449629587435035051,0.045503704125649649492 , 0.007091003462846911073 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Tetrahedron Gauss-Legendre quadrature 4 ";
        return buffer.str();
    }

protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TetrahedronGaussLegendreIntegrationPoints4

class KRATOS_API(KRATOS_CORE) TetrahedronGaussLegendreIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedronGaussLegendreIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 3;

    typedef IntegrationPoint<3> IntegrationPointType;

    typedef std::array<IntegrationPointType, 24> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 24;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.35619138622025439121,0.21460287125991520293,0.21460287125991520293 , 0.006653791709644939366 );
        msIntegrationPoints[1] = IntegrationPointType( 0.21460287125991520293,0.35619138622025439121,0.21460287125991520293 , 0.006653791709644939366 );
        msIntegrationPoints[2] = IntegrationPointType( 0.21460287125991520293,0.21460287125991520293,0.35619138622025439121 , 0.006653791709644939366);
        msIntegrationPoints[3] = IntegrationPointType( 0.21460287125991520293,0.21460287125991520293,0.21460287125991520293 , 0.006653791709644939366 );
        msIntegrationPoints[4] = IntegrationPointType( 0.87797812439616594065,0.040673958534611353116,0.040673958534611353116 , 0.001679535175890970435 );
        msIntegrationPoints[5] = IntegrationPointType( 0.040673958534611353116,0.87797812439616594065,0.040673958534611353116 , 0.001679535175890970435 );
        msIntegrationPoints[6] = IntegrationPointType( 0.040673958534611353116,0.040673958534611353116,0.87797812439616594065 , 0.001679535175890970435 );
        msIntegrationPoints[7] = IntegrationPointType( 0.040673958534611353116,0.040673958534611353116,0.040673958534611353116 , 0.001679535175890970435 );
        msIntegrationPoints[8] = IntegrationPointType( 0.032986329573173468968,0.32233789014227551034,0.32233789014227551034 , 0.009226196923987899723 );
        msIntegrationPoints[9] = IntegrationPointType( 0.32233789014227551034,0.032986329573173468968,0.32233789014227551034 , 0.009226196923987899723 );
        msIntegrationPoints[10] = IntegrationPointType( 0.32233789014227551034,0.32233789014227551034,0.032986329573173468968 , 0.009226196923987899723 );
        msIntegrationPoints[11] = IntegrationPointType( 0.32233789014227551034,0.32233789014227551034,0.32233789014227551034 , 0.009226196923987899723 );
        msIntegrationPoints[12] = IntegrationPointType( 0.60300566479164914137,0.26967233145831580803,0.063661001875017525299 , 0.008035714285714285714 );
        msIntegrationPoints[13] = IntegrationPointType( 0.60300566479164914137,0.063661001875017525299,0.26967233145831580803 , 0.008035714285714285714 );
        msIntegrationPoints[14] = IntegrationPointType( 0.60300566479164914137,0.063661001875017525299,0.063661001875017525299 , 0.008035714285714285714 );
        msIntegrationPoints[15] = IntegrationPointType( 0.063661001875017525299,0.60300566479164914137,0.26967233145831580803 , 0.008035714285714285714 );
        msIntegrationPoints[16] = IntegrationPointType( 0.063661001875017525299,0.60300566479164914137,0.063661001875017525299 , 0.008035714285714285714 );
        msIntegrationPoints[17] = IntegrationPointType( 0.063661001875017525299,0.063661001875017525299,0.60300566479164914137 ,0.008035714285714285714 );
        msIntegrationPoints[18] = IntegrationPointType( 0.26967233145831580803,0.60300566479164914137,0.063661001875017525299 , 0.008035714285714285714 );
        msIntegrationPoints[19] = IntegrationPointType( 0.26967233145831580803,0.063661001875017525299,0.60300566479164914137 , 0.008035714285714285714 );
        msIntegrationPoints[20] = IntegrationPointType( 0.26967233145831580803,0.063661001875017525299,0.063661001875017525299 , 0.008035714285714285714 );
        msIntegrationPoints[21] = IntegrationPointType( 0.063661001875017525299,0.26967233145831580803,0.60300566479164914137 , 0.008035714285714285714 );
        msIntegrationPoints[22] = IntegrationPointType( 0.063661001875017525299,0.26967233145831580803,0.063661001875017525299 , 0.008035714285714285714 );
        msIntegrationPoints[23] = IntegrationPointType( 0.063661001875017525299,0.063661001875017525299,0.26967233145831580803 , 0.008035714285714285714 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Tetrahedron Gauss-Legendre quadrature 5 ";
        return buffer.str();
    }

protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TetrahedronGaussLegendreIntegrationPoints5


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined


