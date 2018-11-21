#include "catch.hpp"

#include <ANurbs/Core>

TEST_CASE( "Spatial B-Spline curve on surface", "[CurveOnSurface3d]" ) {
    using namespace ANurbs;

    Pointer<CurveGeometry2D> curveGeometry;
    {
        int degree = 2;
        int nbPoles = 3;
        bool isRational = false;

        curveGeometry = New<CurveGeometry2D>(degree, nbPoles, isRational);

        curveGeometry->SetKnot( 0, 7.0);
        curveGeometry->SetKnot( 1, 7.0);
        curveGeometry->SetKnot( 2, 9.0);
        curveGeometry->SetKnot( 3, 9.0);

        curveGeometry->SetPole(0, {3, 2});
        curveGeometry->SetPole(1, {1, 4});
        curveGeometry->SetPole(2, {2, 5});
    }

    Pointer<SurfaceGeometry3D> surfaceGeometry;
    {
        int degreeU = 2;
        int degreeV = 2;

        int nbPolesU = 3;
        int nbPolesV = 3;

        bool isRational = false;

        auto knotsU = {1, 1, 3, 3};

        surfaceGeometry = New<SurfaceGeometry3D>(degreeU, degreeV, nbPolesU,
            nbPolesV, isRational);

        surfaceGeometry->SetKnotU(0, 1);
        surfaceGeometry->SetKnotU(1, 1);
        surfaceGeometry->SetKnotU(2, 3);
        surfaceGeometry->SetKnotU(3, 3);

        surfaceGeometry->SetKnotV(0, 2);
        surfaceGeometry->SetKnotV(1, 2);
        surfaceGeometry->SetKnotV(2, 6);
        surfaceGeometry->SetKnotV(3, 6);

        surfaceGeometry->SetPole(0, 0, { 0,  0,  3});
        surfaceGeometry->SetPole(0, 1, { 0,  5,  0});
        surfaceGeometry->SetPole(0, 2, { 0, 10,  2});
        surfaceGeometry->SetPole(1, 0, { 5,  0,  5});
        surfaceGeometry->SetPole(1, 1, { 5,  5,  0});
        surfaceGeometry->SetPole(1, 2, { 5, 10,  3});
        surfaceGeometry->SetPole(2, 0, {10,  0,  1});
        surfaceGeometry->SetPole(2, 1, {10,  5, -1});
        surfaceGeometry->SetPole(2, 2, {10, 10,  0});
    }

    CurveOnSurface3D curveOnSurface(curveGeometry, surfaceGeometry,
        curveGeometry->Domain());

    SECTION( "Check point evaluation" ) {
        auto point = curveOnSurface.PointAt(8);

        REQUIRE( point[0] == Approx(3.75        ) );
        REQUIRE( point[1] == Approx(4.375       ) );
        REQUIRE( point[2] == Approx(1.5063476563) );
    }

    SECTION( "Check derivative evaluation" ) {
        auto point = curveOnSurface.DerivativesAt(8, 2);

        REQUIRE( point[0][0] == Approx( 3.75        ) );
        REQUIRE( point[0][1] == Approx( 4.375       ) );
        REQUIRE( point[0][2] == Approx( 1.5063476563) );

        REQUIRE( point[1][0] == Approx(-2.5         ) );
        REQUIRE( point[1][1] == Approx( 3.75        ) );
        REQUIRE( point[1][2] == Approx(-0.658203125 ) );
    }

    SECTION( "Check spans" ) {
        auto spans = curveOnSurface.Spans();

        REQUIRE( spans.size() == 1 );

        REQUIRE( spans[0].T0() == Approx(7.0) );
        REQUIRE( spans[0].T1() == Approx(9.0) );
    }
}
