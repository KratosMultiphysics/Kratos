#include "catch.hpp"

#include <ANurbs/src/Curve.h>

using namespace ANurbs;


TEST_CASE( "Spatial B-Spline curve", "[Curve]" ) {
    int degree = 4;
    int nbPoles = 8;
    bool isRational = false;

    Pointer<CurveGeometry3D> geometry = New<CurveGeometry3D>(degree, nbPoles,
        isRational);

    geometry->SetKnot( 0,  3.0);
    geometry->SetKnot( 1,  3.0);
    geometry->SetKnot( 2,  3.0);
    geometry->SetKnot( 3,  3.0);
    geometry->SetKnot( 4,  6.5);
    geometry->SetKnot( 5, 10.0);
    geometry->SetKnot( 6, 13.5);
    geometry->SetKnot( 7, 17.0);
    geometry->SetKnot( 8, 17.0);
    geometry->SetKnot( 9, 17.0);
    geometry->SetKnot(10, 17.0);

    geometry->SetPole(0, {  0, -25, - 5});
    geometry->SetPole(1, {-15, -15,   0});
    geometry->SetPole(2, {  5, - 5, - 3});
    geometry->SetPole(3, { 15, -15,   3});
    geometry->SetPole(4, { 25,   0,   6});
    geometry->SetPole(5, { 15,  15,   6});
    geometry->SetPole(6, {- 5, - 5, - 3});
    geometry->SetPole(7, {-25,  15,   4});

    Curve3D curve(geometry);

    SECTION( "Compute arc length" ) {
        double length = curve.Length();

        CHECK( length == Approx(104.186832011) );
    }

    SECTION( "Compute arc length of segment" ) {
        double t0 =  5.0;
        double t1 = 10.0;

        double length = curve.SegmentLength(t0, t1);

        CHECK( length == Approx(24.4282553458) );
    }
}
