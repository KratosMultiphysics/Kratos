#include "catch.hpp"

#include <ANurbs/src/SurfaceGeometry.h>

using namespace ANurbs;


TEST_CASE( "Geometry of a spatial B-Spline surface", "[SurfaceGeometry]" ) {
    int degreeU = 2;
    int degreeV = 1;

    int nbPolesU = 4;
    int nbPolesV = 3;

    bool isRational = false;

    SurfaceGeometry3D geometry(degreeU, degreeV, nbPolesU, nbPolesV,
        isRational);

    geometry.SetKnotU(0,  0.0);
    geometry.SetKnotU(1,  0.0);
    geometry.SetKnotU(2,  7.5);
    geometry.SetKnotU(3, 15.0);
    geometry.SetKnotU(4, 15.0);

    geometry.SetKnotV(0,  0.0);
    geometry.SetKnotV(1, 10.0);
    geometry.SetKnotV(2, 20.0);

    geometry.SetPole(0, 0, {-10.0, - 5.0, -1.0         });
    geometry.SetPole(0, 1, {-12.0,   3.0,  3.0         });
    geometry.SetPole(0, 2, {- 9.0,  11.0, -0.0701928417});
    geometry.SetPole(1, 0, {- 5.0, - 3.0,  1.0         });
    geometry.SetPole(1, 1, {- 6.0,   4.0, -2.0         });
    geometry.SetPole(1, 2, {- 5.0,   7.0,  0.9298071583});
    geometry.SetPole(2, 0, {  0.0, - 4.0, -1.0         });
    geometry.SetPole(2, 1, {  1.0,   6.0,  5.0         });
    geometry.SetPole(2, 2, {  0.0,  13.0, -0.2350184214});
    geometry.SetPole(3, 0, {  4.0, - 2.0,  0.0         });
    geometry.SetPole(3, 1, {  5.0,   4.0, -1.0         });
    geometry.SetPole(3, 2, {  5.0,  11.0,  0.7649815786});

    SECTION( "Check Poles" ) {
        const auto poles = geometry.Poles();

        for (int i = 0; i < poles.NbRows(); i++) {
            for (int j = 0; j < poles.NbCols(); j++) {
                const auto actual = poles(i, j);
                const auto expected = geometry.Pole(i, j);
                CHECK( actual[0] == expected[0] );
                CHECK( actual[1] == expected[1] );
                CHECK( actual[2] == expected[2] );
            }
        }
    }

    SECTION( "Check Weights" ) {
        const auto weights = geometry.Weights();

        for (int i = 0; i < weights.NbRows(); i++) {
            for (int j = 0; j < weights.NbCols(); j++) {
                const auto actual = weights(i, j);
                const auto expected = geometry.Weight(i, j);
                CHECK( actual == expected );
            }
        }
    }

    SECTION( "Degrees" ) {
        REQUIRE( geometry.DegreeU() == 2 );
        REQUIRE( geometry.DegreeV() == 1 );
    }

    SECTION( "Number of knots" ) {
        REQUIRE( geometry.NbKnotsU() == 5 );
        REQUIRE( geometry.NbKnotsV() == 3 );
    }

    SECTION( "Number of poles" ) {
        REQUIRE( geometry.NbPolesU() ==  4 );
        REQUIRE( geometry.NbPolesV() ==  3 );
        REQUIRE( geometry.NbPoles()  == 12 );
    }

    SECTION( "Evaluate point" ) {
        double u = 12;
        double v = 5;

        auto point = geometry.PointAt(u, v);

        REQUIRE( point[0] == Approx( 1.46) );
        REQUIRE( point[1] == Approx( 0.96) );
        REQUIRE( point[2] == Approx( 0.9 ) );
    }

    SECTION( "Evaluate derivatives" ) {
        double u = 12;
        double v = 5;

        int order = 2;

        auto derivatives = geometry.DerivativesAt(u, v, order);

        REQUIRE( derivatives.size() == 6 );

        REQUIRE( derivatives[0][0] == Approx( 1.46        ) );
        REQUIRE( derivatives[0][1] == Approx( 0.96        ) );
        REQUIRE( derivatives[0][2] == Approx( 0.9         ) );

        REQUIRE( derivatives[1][0] == Approx( 0.96        ) );
        REQUIRE( derivatives[1][1] == Approx( 0.0266666667) );
        REQUIRE( derivatives[1][2] == Approx(-0.2666666667) );

        REQUIRE( derivatives[2][0] == Approx( 0.084       ) );
        REQUIRE( derivatives[2][1] == Approx( 0.832       ) );
        REQUIRE( derivatives[2][2] == Approx( 0.276       ) );

        REQUIRE( derivatives[3][0] == Approx( 0.0355555556) );
        REQUIRE( derivatives[3][1] == Approx(-0.0088888889) );
        REQUIRE( derivatives[3][2] == Approx(-0.1333333333) );

        REQUIRE( derivatives[4][0] == Approx( 0.0106666667) );
        REQUIRE( derivatives[4][1] == Approx(-0.048       ) );
        REQUIRE( derivatives[4][2] == Approx(-0.064       ) );

        REQUIRE( derivatives[5][0] == Approx( 0           ) );
        REQUIRE( derivatives[5][1] == Approx( 0           ) );
        REQUIRE( derivatives[5][2] == Approx( 0           ) );
    }
}

TEST_CASE( "Geometry of a spatial Nurbs surface", "[SurfaceGeometry]" ) {
    int degreeU = 2;
    int degreeV = 1;

    int nbPolesU = 4;
    int nbPolesV = 3;

    bool isRational = true;

    SurfaceGeometry3D geometry(degreeU, degreeV, nbPolesU, nbPolesV,
        isRational);

    geometry.SetKnotU(0,  0.0);
    geometry.SetKnotU(1,  0.0);
    geometry.SetKnotU(2,  7.5);
    geometry.SetKnotU(3, 15.0);
    geometry.SetKnotU(4, 15.0);

    geometry.SetKnotV(0,  0.0);
    geometry.SetKnotV(1, 10.0);
    geometry.SetKnotV(2, 20.0);

    geometry.SetPole(0, 0, {-10.0, - 5.0, -1.0         });
    geometry.SetPole(0, 1, {-12.0,   3.0,  3.0         });
    geometry.SetPole(0, 2, {- 9.0,  11.0, -0.0701928417});
    geometry.SetPole(1, 0, {- 5.0, - 3.0,  1.0         });
    geometry.SetPole(1, 1, {- 6.0,   4.0, -2.0         });
    geometry.SetPole(1, 2, {- 5.0,   7.0,  0.9298071583});
    geometry.SetPole(2, 0, {  0.0, - 4.0, -1.0         });
    geometry.SetPole(2, 1, {  1.0,   6.0,  5.0         });
    geometry.SetPole(2, 2, {  0.0,  13.0, -0.2350184214});
    geometry.SetPole(3, 0, {  4.0, - 2.0,  0.0         });
    geometry.SetPole(3, 1, {  5.0,   4.0, -1.0         });
    geometry.SetPole(3, 2, {  5.0,  11.0,  0.7649815786});

    geometry.SetWeight(0, 0, 1.0);
    geometry.SetWeight(0, 1, 1.0);
    geometry.SetWeight(0, 2, 1.0);
    geometry.SetWeight(1, 0, 1.0);
    geometry.SetWeight(1, 1, 2.5);
    geometry.SetWeight(1, 2, 1.0);
    geometry.SetWeight(2, 0, 1.0);
    geometry.SetWeight(2, 1, 1.0);
    geometry.SetWeight(2, 2, 1.0);
    geometry.SetWeight(3, 0, 1.0);
    geometry.SetWeight(3, 1, 1.0);
    geometry.SetWeight(3, 2, 1.0);

    SECTION( "Degrees" ) {
        REQUIRE( geometry.DegreeU() == 2 );
        REQUIRE( geometry.DegreeV() == 1 );
    }

    SECTION( "Number of knots" ) {
        REQUIRE( geometry.NbKnotsU() == 5 );
        REQUIRE( geometry.NbKnotsV() == 3 );
    }

    SECTION( "Number of poles" ) {
        REQUIRE( geometry.NbPolesU() ==  4 );
        REQUIRE( geometry.NbPolesV() ==  3 );
        REQUIRE( geometry.NbPoles()  == 12 );
    }

    SECTION( "Evaluate point" ) {
        double u = 12;
        double v = 5;

        auto point = geometry.PointAt(u, v);

        REQUIRE( point[0] == Approx( 1.0377358491) );
        REQUIRE( point[1] == Approx( 1.1320754717) );
        REQUIRE( point[2] == Approx( 0.7358490566) );
    }

    SECTION( "Evaluate derivatives" ) {
        double u = 12;
        double v = 5;

        int order = 2;

        auto derivatives = geometry.DerivativesAt(u, v, order);

        REQUIRE( derivatives.size() == 6 );

        REQUIRE( derivatives[0][0] == Approx( 1.0377358491) );
        REQUIRE( derivatives[0][1] == Approx( 1.1320754717) );
        REQUIRE( derivatives[0][2] == Approx( 0.7358490566) );

        REQUIRE( derivatives[1][0] == Approx( 1.1712353151) );
        REQUIRE( derivatives[1][1] == Approx(-0.0830663344) );
        REQUIRE( derivatives[1][2] == Approx(-0.14833274  ) );

        REQUIRE( derivatives[2][0] == Approx(-0.0004271983) );
        REQUIRE( derivatives[2][1] == Approx( 0.8173727305) );
        REQUIRE( derivatives[2][2] == Approx( 0.2294054824) );

        REQUIRE( derivatives[3][0] == Approx( 0.0334131158) );
        REQUIRE( derivatives[3][1] == Approx( 0.0214196208) );
        REQUIRE( derivatives[3][2] == Approx(-0.171394283 ) );

        REQUIRE( derivatives[4][0] == Approx( 0.0499024922) );
        REQUIRE( derivatives[4][1] == Approx(-0.0351431047) );
        REQUIRE( derivatives[4][2] == Approx(-0.0293933919) );

        REQUIRE( derivatives[5][0] == Approx( 9.6724e-6   ) );
        REQUIRE( derivatives[5][1] == Approx(-0.0185065524) );
        REQUIRE( derivatives[5][2] == Approx(-0.0051940864) );
    }
}
