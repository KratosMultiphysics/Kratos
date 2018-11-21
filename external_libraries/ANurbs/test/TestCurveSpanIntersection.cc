#include "catch.hpp"

#include <ANurbs/src/Curve.h>
#include <ANurbs/src/CurveGeometry.h>
#include <ANurbs/src/CurveSpanIntersection.h>

TEST_CASE( "Intersect curve with surface knot spans (1)",
    "[CurveSpanIntersection][Curve][CurveGeometry][2D]" ) {
    using namespace ANurbs;

    // create curve

    Pointer<CurveGeometry2D> geometry = New<CurveGeometry2D>(
        2,      // degree
        5,      // nbPoles
        false   // isRational
    );

    geometry->SetKnot(0,  1.0         );
    geometry->SetKnot(1,  1.0         );
    geometry->SetKnot(2,  4.3333333333);
    geometry->SetKnot(3,  7.6666666667);
    geometry->SetKnot(4, 11.0         );
    geometry->SetKnot(5, 11.0         );

    geometry->SetPole(0, { 5.0,  5.0});
    geometry->SetPole(1, { 8.0,  7.0});
    geometry->SetPole(2, { 6.0,  8.0});
    geometry->SetPole(3, { 8.0, 10.0});
    geometry->SetPole(4, { 9.0,  8.0});

    Curve2D curve(geometry);

    // spans

    std::vector<double> knotsU = {5, 5, 6, 7, 8, 9, 10, 10};
    std::vector<double> knotsV = {5, 5, 6, 8, 8, 10, 10};

    // compute intersections

    CurveSpanIntersection2D intersection;
    intersection.Compute(curve, knotsU, knotsV, 1e-7, false);

    SECTION( "Check intersections" ) {
        CHECK( intersection.NbIntersections() == 9 );

        CHECK( intersection.Parameter(0) == Approx( 1.0         ) );
        CHECK( intersection.Parameter(1) == Approx( 1.6366100188) );
        CHECK( intersection.Parameter(2) == Approx( 1.9308025998) );
        CHECK( intersection.Parameter(3) == Approx( 2.6666666667) );
        CHECK( intersection.Parameter(4) == Approx( 4.3333333333) );
        CHECK( intersection.Parameter(5) == Approx( 5.7140452079) );
        CHECK( intersection.Parameter(6) == Approx( 7.6666666667) );
        CHECK( intersection.Parameter(7) == Approx( 9.3333333333) );
        CHECK( intersection.Parameter(8) == Approx(11.0         ) );
    }
}

TEST_CASE( "Intersect curve with surface knot spans (2)",
    "[CurveSpanIntersection][Curve][CurveGeometry][2D]" ) {
    using namespace ANurbs;

    // create curve

    Pointer<CurveGeometry2D> geometry = New<CurveGeometry2D>(
        2,      // degree
        5,      // nbPoles
        false   // isRational
    );

    geometry->SetKnot(0, 2.0);
    geometry->SetKnot(1, 2.0);
    geometry->SetKnot(2, 3.0);
    geometry->SetKnot(3, 3.0);
    geometry->SetKnot(4, 4.0);
    geometry->SetKnot(5, 4.0);

    geometry->SetPole(0, { 8.0,  6.0});
    geometry->SetPole(1, { 6.0,  7.0});
    geometry->SetPole(2, { 8.0,  8.0});
    geometry->SetPole(3, {10.0,  9.0});
    geometry->SetPole(4, { 8.0, 10.0});

    Curve2D curve(geometry);

    // spans

    std::vector<double> knotsU = {5, 5, 6, 7, 8, 9, 10, 10};
    std::vector<double> knotsV = {5, 5, 6, 8, 8, 10, 10};

    // compute intersections

    CurveSpanIntersection2D intersection;
    intersection.Compute(curve, knotsU, knotsV, 1e-7, false);

    SECTION( "Check intersections" ) {
        CHECK( intersection.NbIntersections() == 5 );

        CHECK( intersection.Parameter(0) == Approx(2.0) );
        CHECK( intersection.Parameter(1) == Approx(2.5) );
        CHECK( intersection.Parameter(2) == Approx(3.0) );
        CHECK( intersection.Parameter(3) == Approx(3.5) );
        CHECK( intersection.Parameter(4) == Approx(4.0) );
    }
}

TEST_CASE( "Intersect curve with surface knot spans (3)",
    "[CurveSpanIntersection][Curve][CurveGeometry][2D]" ) {
    using namespace ANurbs;

    // create curve

    Pointer<CurveGeometry2D> geometry = New<CurveGeometry2D>(
        2,      // degree
        4,      // nbPoles
        false   // isRational
    );

    geometry->SetKnot(0, 1.0);
    geometry->SetKnot(1, 1.0);
    geometry->SetKnot(2, 2.0);
    geometry->SetKnot(3, 3.0);
    geometry->SetKnot(4, 3.0);

    geometry->SetPole(0, {9.0, 6.0});
    geometry->SetPole(1, {9.0, 7.0});
    geometry->SetPole(2, {6.0, 7.0});
    geometry->SetPole(3, {6.0, 8.0});

    Curve2D curve(geometry);

    // spans

    std::vector<double> knotsU = {5, 5, 6, 7, 8, 9, 10, 10};
    std::vector<double> knotsV = {5, 5, 6, 8, 8, 10, 10};

    // compute intersections

    CurveSpanIntersection2D intersection;
    intersection.Compute(curve, knotsU, knotsV, 1e-7, false);

    SECTION( "Check intersections" ) {
        CHECK( intersection.NbIntersections() == 4 );

        CHECK( intersection.Parameter(0) == Approx(1.0         ) );
        CHECK( intersection.Parameter(1) == Approx(1.8164965809) );
        CHECK( intersection.Parameter(2) == Approx(2.1835034191) );
        CHECK( intersection.Parameter(3) == Approx(3.0         ) );
    }
}

TEST_CASE( "Intersect curve with surface knot spans (4)",
    "[CurveSpanIntersection][Curve][CurveGeometry][2D]" ) {
    using namespace ANurbs;

    // create curve

    Pointer<CurveGeometry2D> geometry = New<CurveGeometry2D>(
        2,      // degree
        4,      // nbPoles
        false   // isRational
    );

    geometry->SetKnot(0, 1.0);
    geometry->SetKnot(1, 1.0);
    geometry->SetKnot(2, 2.0);
    geometry->SetKnot(3, 3.0);
    geometry->SetKnot(4, 3.0);

    geometry->SetPole(0, {9.0, 6.0});
    geometry->SetPole(1, {9.0, 7.0});
    geometry->SetPole(2, {6.0, 7.0});
    geometry->SetPole(3, {6.0, 8.0});

    Curve2D curve(geometry);

    // spans

    std::vector<double> knotsU = {5, 5, 6, 7, 8, 9, 10, 10};
    std::vector<double> knotsV = {5, 5, 6, 8, 8, 10, 10};

    // compute intersections

    CurveSpanIntersection2D intersection;
    intersection.Compute(curve, knotsU, knotsV, 1e-7, true);

    SECTION( "Check intersections" ) {
        CHECK( intersection.NbIntersections() == 5 );

        CHECK( intersection.Parameter(0) == Approx(1.0         ) );
        CHECK( intersection.Parameter(1) == Approx(1.8164965809) );
        CHECK( intersection.Parameter(2) == Approx(2.0         ) );
        CHECK( intersection.Parameter(3) == Approx(2.1835034191) );
        CHECK( intersection.Parameter(4) == Approx(3.0         ) );
    }
}
