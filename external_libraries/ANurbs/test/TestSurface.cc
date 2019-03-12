#include "catch.hpp"

#include <ANurbs/src/Surface.h>

using namespace ANurbs;


TEST_CASE( "Geometry of a spatial B-Spline srf", "[Surface]" ) {
    int degreeU = 2;
    int degreeV = 1;

    int nbPolesU = 4;
    int nbPolesV = 3;

    bool isRational = false;

    Pointer<SurfaceGeometry3D> geometry = New<SurfaceGeometry3D>(degreeU,
        degreeV, nbPolesU, nbPolesV, isRational);

    geometry->SetKnotU(0,  0.0);
    geometry->SetKnotU(1,  0.0);
    geometry->SetKnotU(2,  7.5);
    geometry->SetKnotU(3, 15.0);
    geometry->SetKnotU(4, 15.0);

    geometry->SetKnotV(0,  0.0);
    geometry->SetKnotV(1, 10.0);
    geometry->SetKnotV(2, 20.0);

    geometry->SetPole(0, 0, {-10.0, - 5.0, -1.0         });
    geometry->SetPole(0, 1, {-12.0,   3.0,  3.0         });
    geometry->SetPole(0, 2, {- 9.0,  11.0, -0.0701928417});
    geometry->SetPole(1, 0, {- 5.0, - 3.0,  1.0         });
    geometry->SetPole(1, 1, {- 6.0,   4.0, -2.0         });
    geometry->SetPole(1, 2, {- 5.0,   7.0,  0.9298071583});
    geometry->SetPole(2, 0, {  0.0, - 4.0, -1.0         });
    geometry->SetPole(2, 1, {  1.0,   6.0,  5.0         });
    geometry->SetPole(2, 2, {  0.0,  13.0, -0.2350184214});
    geometry->SetPole(3, 0, {  4.0, - 2.0,  0.0         });
    geometry->SetPole(3, 1, {  5.0,   4.0, -1.0         });
    geometry->SetPole(3, 2, {  6.0,  11.0,  0.7649815786});

    Surface3D surface(geometry);

    SECTION( "Evaluate point" ) {
        double u = 12;
        double v = 5;

        auto point = geometry->PointAt(u, v);

        REQUIRE( point[0] == Approx( 1.46) );
        REQUIRE( point[1] == Approx( 0.96) );
        REQUIRE( point[2] == Approx( 0.9 ) );
    }

    SECTION( "Evaluate derivatives" ) {
        double u = 12;
        double v = 5;

        int order = 2;

        auto derivatives = geometry->DerivativesAt(u, v, order);

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
