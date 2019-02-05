//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    hbui
//
//


#if !defined(KRATOS_TRIANGLE_GAUSS_RADAU_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_TRIANGLE_GAUSS_RADAU_INTEGRATION_POINTS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "integration/quadrature.h"

/// REF: G. R. Cowper, GAUSSIAN QUADRATURE FORMULAS FOR TRIANGLES

namespace Kratos
{

class KRATOS_API(KRATOS_CORE) TriangleGaussRadauIntegrationPoints1
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussRadauIntegrationPoints1);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 3> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.5 , 0.5 , 1.00 / 3.00 );
        msIntegrationPoints[1] = IntegrationPointType( 0.5 , 0.0 , 1.00 / 3.00 );
        msIntegrationPoints[2] = IntegrationPointType( 0.0 , 0.5 , 1.00 / 3.00 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Radau quadrature 1 (3 points, degree of precision = 2) ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussRadauIntegrationPoints1


class KRATOS_API(KRATOS_CORE) TriangleGaussRadauIntegrationPoints2
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussRadauIntegrationPoints2);
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
        msIntegrationPoints[0] = IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -0.562500000000000 );
        msIntegrationPoints[1] = IntegrationPointType( 0.6 , 0.2 , 0.520833333333333 );
        msIntegrationPoints[2] = IntegrationPointType( 0.2 , 0.6 , 0.520833333333333 );
        msIntegrationPoints[3] = IntegrationPointType( 0.2 , 0.2 , 0.520833333333333 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Radau quadrature 2 (4 points, degree of precision = 3)";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussRadauIntegrationPoints2


class KRATOS_API(KRATOS_CORE) TriangleGaussRadauIntegrationPoints3
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussRadauIntegrationPoints3);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 6> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.816847572980459 , 0.091576213509771 , 0.109951743655322 );
        msIntegrationPoints[1] = IntegrationPointType( 0.091576213509771 , 0.816847572980459 , 0.109951743655322 );
        msIntegrationPoints[2] = IntegrationPointType( 0.091576213509771 , 0.091576213509771 , 0.109951743655322 );
        msIntegrationPoints[3] = IntegrationPointType( 0.108103018168070 , 0.445948490915965 , 0.223381589678011 );
        msIntegrationPoints[4] = IntegrationPointType( 0.445948490915965 , 0.108103018168070 , 0.223381589678011 );
        msIntegrationPoints[5] = IntegrationPointType( 0.445948490915965 , 0.445948490915965 , 0.223381589678011 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Radau quadrature 3 (6 points, degree of precision = 4) ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussRadauIntegrationPoints2


class KRATOS_API(KRATOS_CORE) TriangleGaussRadauIntegrationPoints4
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussRadauIntegrationPoints4);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 7> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 7;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , 0.225000000000000 );
        msIntegrationPoints[1] = IntegrationPointType( 0.797426985353087 , 0.101286507323456 , 0.125939180544827 );
        msIntegrationPoints[2] = IntegrationPointType( 0.101286507323456 , 0.797426985353087 , 0.125939180544827 );
        msIntegrationPoints[3] = IntegrationPointType( 0.101286507323456 , 0.101286507323456 , 0.125939180544827 );
        msIntegrationPoints[4] = IntegrationPointType( 0.470142064105115 , 0.059715871789770 , 0.132394152788506 );
        msIntegrationPoints[5] = IntegrationPointType( 0.059715871789770 , 0.470142064105115 , 0.132394152788506 );
        msIntegrationPoints[6] = IntegrationPointType( 0.470142064105115 , 0.470142064105115 , 0.132394152788506 );
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Radau quadrature 4 (7 points, degree of precision = 5) ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussRadauIntegrationPoints4


class KRATOS_API(KRATOS_CORE) TriangleGaussRadauIntegrationPoints5
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussRadauIntegrationPoints5);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 12> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 12;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 0.873821971016996 , 0.063089014491502 , 0.050844906370207 );
        msIntegrationPoints[1] = IntegrationPointType( 0.063089014491502 , 0.873821971016996 , 0.050844906370207 );
        msIntegrationPoints[2] = IntegrationPointType( 0.063089014491502 , 0.063089014491502 , 0.050844906370207 );

        msIntegrationPoints[3] = IntegrationPointType( 0.501426509658179 , 0.249286745170910 , 0.116786275726379 );
        msIntegrationPoints[4] = IntegrationPointType( 0.249286745170910 , 0.501426509658179 , 0.116786275726379 );
        msIntegrationPoints[5] = IntegrationPointType( 0.249286745170910 , 0.249286745170910 , 0.116786275726379 );

        msIntegrationPoints[6] = IntegrationPointType( 0.636502499121399 , 0.310352451033785 , 0.082851075618374 );
        msIntegrationPoints[7] = IntegrationPointType( 0.636502499121399 , 0.053145049844816 , 0.082851075618374 );
        msIntegrationPoints[8] = IntegrationPointType( 0.310352451033785 , 0.636502499121399 , 0.082851075618374 );
        msIntegrationPoints[9] = IntegrationPointType( 0.310352451033785 , 0.053145049844816 , 0.082851075618374 );
        msIntegrationPoints[10] = IntegrationPointType( 0.053145049844816 , 0.636502499121399 , 0.082851075618374 );
        msIntegrationPoints[11] = IntegrationPointType( 0.053145049844816 , 0.310352451033785 , 0.082851075618374 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Radau quadrature 5 (12 points, degree of precision = 6) ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussRadauIntegrationPoints5


class KRATOS_API(KRATOS_CORE) TriangleGaussRadauIntegrationPoints6
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TriangleGaussRadauIntegrationPoints6);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef std::array<IntegrationPointType, 13> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return 13;
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints[0] = IntegrationPointType( 1.00 / 3.00 , 1.00 / 3.00 , -0.149570044467670 );

        msIntegrationPoints[1] = IntegrationPointType( 0.479308067841923 , 0.260345966079038 , 0.175615257433204 );
        msIntegrationPoints[2] = IntegrationPointType( 0.260345966079038 , 0.479308067841923 , 0.175615257433204 );
        msIntegrationPoints[3] = IntegrationPointType( 0.260345966079038 , 0.260345966079038 , 0.175615257433204 );

        msIntegrationPoints[4] = IntegrationPointType( 0.869739794195568 , 0.065130102902216 , 0.053347235608839 );
        msIntegrationPoints[5] = IntegrationPointType( 0.065130102902216 , 0.869739794195568 , 0.053347235608839 );
        msIntegrationPoints[6] = IntegrationPointType( 0.065130102902216 , 0.065130102902216 , 0.053347235608839 );

        msIntegrationPoints[7] = IntegrationPointType( 0.638444188569809 , 0.312865496004875 , 0.077113760890257 );
        msIntegrationPoints[8] = IntegrationPointType( 0.638444188569809 , 0.048690315425316 , 0.077113760890257 );
        msIntegrationPoints[9] = IntegrationPointType( 0.312865496004875 , 0.638444188569809 , 0.077113760890257 );
        msIntegrationPoints[10] = IntegrationPointType( 0.312865496004875 , 0.048690315425316 , 0.077113760890257 );
        msIntegrationPoints[11] = IntegrationPointType( 0.048690315425316 , 0.638444188569809 , 0.077113760890257 );
        msIntegrationPoints[12] = IntegrationPointType( 0.048690315425316 , 0.312865496004875 , 0.077113760890257 );

        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Triangle Gauss-Radau quadrature 6 (13 points, degree of precision = 7) ";
        return buffer.str();
    }
protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;

}; // Class TriangleGaussRadauIntegrationPoints6


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_TRIANGLE_GAUSS_RADAU_INTEGRATION_POINTS_H_INCLUDED  defined


