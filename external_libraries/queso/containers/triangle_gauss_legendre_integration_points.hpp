// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H
#define TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H

//// STL includes
#include <vector>
//// Projecy includes
#include "containers/integration_point.hpp"

namespace queso {

class TriangleGaussLegendrePoints1
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 1;
    }

    static const PointArrayType& IntegrationPoints(){

        static const PointArrayType s_integration_points{{
            IntegrationPoint( 1.00 / 3.00 , 1.00 / 3.00 , 1.00 / 2.00 )
        }};

        return s_integration_points;
    }
};

class TriangleGaussLegendrePoints2
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 3;
    }

    static const PointArrayType& IntegrationPoints(){

        static const PointArrayType s_integration_points{{
            IntegrationPoint( 1.00 / 6.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
            IntegrationPoint( 2.00 / 3.00 , 1.00 / 6.00 , 1.00 / 6.00 ),
            IntegrationPoint( 1.00 / 6.00 , 2.00 / 3.00 , 1.00 / 6.00 )
        }};

        return s_integration_points;
    }
};

class TriangleGaussLegendrePoints3
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 6;
    }

    static const PointArrayType& IntegrationPoints(){
        const double wa = 0.054975871827661;
        const double wb = 0.1116907948390055;
        const double Na1 = 0.816847572980459;
        const double Nb1 = 0.108103018168070;
        const double Na2 = 0.091576213509771;
        const double Nb2 = 0.445948490915965;

        static const PointArrayType s_integration_points{{
            IntegrationPoint( Na2, Na2, wa ),
            IntegrationPoint( Na1, Na2, wa ),
            IntegrationPoint( Na2, Na1, wa ),
            IntegrationPoint( Nb2, Nb2, wb ),
            IntegrationPoint( Nb1, Nb2, wb ),
            IntegrationPoint( Nb2, Nb1, wb )
        }};

        return s_integration_points;
    }
};

class TriangleGaussLegendrePoints4
{
public:
    typedef std::size_t SizeType;
    typedef std::vector<IntegrationPoint> PointArrayType;

    static SizeType IntegrationPointsNumber()
    {
        return 12;
    }

    static const PointArrayType& IntegrationPoints(){
        const double wa = 0.025422453185103408460;
        const double wb = 0.058393137863189683013;
        const double wc = 0.041425537809186787597;

        const double N1 = 0.87382197101699554332;
        const double N2 = 0.063089014491502228340;
        const double N3 = 0.50142650965817915742;
        const double N4 = 0.24928674517091042129;
        const double N5 = 0.053145049844816947353;
        const double N6 = 0.31035245103378440542;
        const double N7 = 0.63650249912139864723;

        static const PointArrayType s_integration_points{{
            IntegrationPoint( N1, N2, wa ),
            IntegrationPoint( N2, N1, wa ),
            IntegrationPoint( N2, N2, wa ),
            IntegrationPoint( N3, N4, wb ),
            IntegrationPoint( N4, N3, wb ),
            IntegrationPoint( N4, N4, wb ),
            IntegrationPoint( N5, N6, wc ),
            IntegrationPoint( N6, N5, wc ),
            IntegrationPoint( N5, N7, wc ),
            IntegrationPoint( N6, N7, wc ),
            IntegrationPoint( N7, N5, wc ),
            IntegrationPoint( N7, N6, wc )
        }};
        return s_integration_points;
    }
};
} // End namespace queso

#endif // TRIANGLE_GAUSS_LEGENDRE_INTEGRATION_POINTS_H