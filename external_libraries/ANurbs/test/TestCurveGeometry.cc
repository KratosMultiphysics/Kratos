#include "catch.hpp"

#include <ANurbs/src/CurveGeometry.h>

using namespace ANurbs;


TEST_CASE( "Geometry of a spatial B-Spline curve", "[CurveGeometry]" ) {
    int degree = 4;
    int nbPoles = 8;
    bool I = false;

    CurveGeometry3D geometry(degree, nbPoles, I);

    geometry.SetKnot( 0,   0              );
    geometry.SetKnot( 1,   0              );
    geometry.SetKnot( 2,   0              );
    geometry.SetKnot( 3,   0              );
    geometry.SetKnot( 4,  32.9731425998736);
    geometry.SetKnot( 5,  65.9462851997473);
    geometry.SetKnot( 6,  98.9194277996209);
    geometry.SetKnot( 7, 131.892570399495 );
    geometry.SetKnot( 8, 131.892570399495 );
    geometry.SetKnot( 9, 131.892570399495 );
    geometry.SetKnot(10, 131.892570399495 );

    geometry.SetPole(0, {  0, -25, - 5});
    geometry.SetPole(1, {-15, -15,   0});
    geometry.SetPole(2, {  5, - 5, - 3});
    geometry.SetPole(3, { 15, -15,   3});
    geometry.SetPole(4, { 25,   0,   6});
    geometry.SetPole(5, { 15,  15,   6});
    geometry.SetPole(6, {- 5, - 5, - 3});
    geometry.SetPole(7, {-25,  15,   4});

    REQUIRE( geometry.Dimension()  == 3     );
    REQUIRE( geometry.IsRational() == false );

    REQUIRE( geometry.Degree()  ==  4 );
    REQUIRE( geometry.NbKnots() == 11 );
    REQUIRE( geometry.NbPoles() ==  8 );

    REQUIRE( geometry.Domain().T0() ==   0              );
    REQUIRE( geometry.Domain().T1() == 131.892570399495 );

    SECTION( "Check Poles" ) {
        const auto poles = geometry.Poles();

        for (int i = 0; i < poles.size(); i++) {
            const auto actual = poles[i];
            const auto expected = geometry.Pole(i);
            CHECK( actual[0] == expected[0] );
            CHECK( actual[1] == expected[1] );
            CHECK( actual[2] == expected[2] );
        }
    }

    SECTION( "Check Weights" ) {
        const auto weights = geometry.Weights();

        for (int i = 0; i < weights.size(); i++) {
            const auto actual = weights[i];
            const auto expected = geometry.Weight(i);
            CHECK( actual == expected );
        }
    }

    SECTION( "PointAt(t=0.0)" ) {
        auto point = geometry.PointAt(0.0);

        CHECK( point[0] == Approx(  0) );
        CHECK( point[1] == Approx(-25) );
        CHECK( point[2] == Approx(- 5) );
    }

    SECTION( "PointAt(t=65.9462851997)" ) {
        auto point = geometry.PointAt(65.9462851997);

        CHECK( point[0] == Approx( 18.8888888889) );
        CHECK( point[1] == Approx(- 6.1111111111) );
        CHECK( point[2] == Approx(  4.1666666667) );
    }

    SECTION( "PointAt(t=131.892570399495)" ) {
        auto point = geometry.PointAt(131.892570399495);

        CHECK( point[0] == Approx(-25) );
        CHECK( point[1] == Approx( 15) );
        CHECK( point[2] == Approx(  4) );
    }

    SECTION( "DerivativesAt(t=0.0)" ) {
        auto derivatives = geometry.DerivativesAt(0.0, 3);

        CHECK( derivatives[0][0] == Approx(  0) );
        CHECK( derivatives[0][1] == Approx(-25) );
        CHECK( derivatives[0][2] == Approx(- 5) );

        CHECK( derivatives[1][0] == Approx(- 1.81966277  ) );
        CHECK( derivatives[1][1] == Approx(  1.2131085134) );
        CHECK( derivatives[1][2] == Approx(  0.6065542567) );

        CHECK( derivatives[2][0] == Approx(  0.2759310497) );
        CHECK( derivatives[2][1] == Approx(- 0.0551862099) );
        CHECK( derivatives[2][2] == Approx(- 0.0717420729) );

        CHECK( derivatives[3][0] == Approx(- 0.0189682773) );
        CHECK( derivatives[3][1] == Approx(  0.0005578905) );
        CHECK( derivatives[3][2] == Approx(  0.005523116 ) );
    }

    SECTION( "DerivativesAt(t=65.9462851997)" ) {
        auto derivatives = geometry.DerivativesAt(65.9462851997, 3);

        CHECK( derivatives[0][0] == Approx( 18.8888888889) );
        CHECK( derivatives[0][1] == Approx(- 6.1111111111) );
        CHECK( derivatives[0][2] == Approx(  4.1666666667) );

        CHECK( derivatives[1][0] == Approx(  0.2021847522) );
        CHECK( derivatives[1][1] == Approx(  0.336974587 ) );
        CHECK( derivatives[1][2] == Approx(  0.1010923761) );

        CHECK( derivatives[2][0] == Approx(- 0.0122636022) );
        CHECK( derivatives[2][1] == Approx(  0.0153295028) );
        CHECK( derivatives[2][2] == Approx(- 0.0036790807) );

        CHECK( derivatives[3][0] == Approx(- 0.0005578905) );
        CHECK( derivatives[3][1] == Approx(- 0.0006508723) );
        CHECK( derivatives[3][2] == Approx(  0.0000557891) );
    }

    SECTION( "DerivativesAt(t=131.892570399495)" ) {
        auto derivatives = geometry.DerivativesAt(131.892570399495, 3);

        CHECK( derivatives[0][0] == Approx(-25           ) );
        CHECK( derivatives[0][1] == Approx( 15           ) );
        CHECK( derivatives[0][2] == Approx(  4           ) );

        CHECK( derivatives[1][0] == Approx(- 2.4262170267) );
        CHECK( derivatives[1][1] == Approx(  2.4262170267) );
        CHECK( derivatives[1][2] == Approx(  0.8491759593) );

        CHECK( derivatives[2][0] == Approx(- 0.1103724199) );
        CHECK( derivatives[2][1] == Approx(  0.3311172597) );
        CHECK( derivatives[2][2] == Approx(  0.1269282829) );

        CHECK( derivatives[3][0] == Approx(- 0.0044631241) );
        CHECK( derivatives[3][1] == Approx(  0.0251050729) );
        CHECK( derivatives[3][2] == Approx(  0.0092051934) );
    }
}

TEST_CASE( "Geometry of a spatial Nurbs curve", "[CurveGeometry]" ) {
    int degree = 4;
    int nbPoles = 8;
    bool isRational = true;

    CurveGeometry3D geometry(degree, nbPoles, isRational);

    geometry.SetKnot( 0,   0              );
    geometry.SetKnot( 1,   0              );
    geometry.SetKnot( 2,   0              );
    geometry.SetKnot( 3,   0              );
    geometry.SetKnot( 4,  32.9731425998736);
    geometry.SetKnot( 5,  65.9462851997473);
    geometry.SetKnot( 6,  98.9194277996209);
    geometry.SetKnot( 7, 131.892570399495 );
    geometry.SetKnot( 8, 131.892570399495 );
    geometry.SetKnot( 9, 131.892570399495 );
    geometry.SetKnot(10, 131.892570399495 );

    geometry.SetPole(0, {  0, -25, - 5});
    geometry.SetPole(1, {-15, -15,   0});
    geometry.SetPole(2, {  5, - 5, - 3});
    geometry.SetPole(3, { 15, -15,   3});
    geometry.SetPole(4, { 25,   0,   6});
    geometry.SetPole(5, { 15,  15,   6});
    geometry.SetPole(6, {- 5, - 5, - 3});
    geometry.SetPole(7, {-25,  15,   4});

    geometry.SetWeight(0, 1.0);
    geometry.SetWeight(1, 1.0);
    geometry.SetWeight(2, 1.0);
    geometry.SetWeight(3, 1.0);
    geometry.SetWeight(4, 1.0);
    geometry.SetWeight(5, 1.0);
    geometry.SetWeight(6, 1.0);
    geometry.SetWeight(7, 1.0);

    CHECK( geometry.Dimension()  == 3     );
    CHECK( geometry.IsRational() == true  );

    CHECK( geometry.Degree()  ==  4 );
    CHECK( geometry.NbKnots() == 11 );
    CHECK( geometry.NbPoles() ==  8 );

    CHECK( geometry.Domain().T0() ==   0              );
    CHECK( geometry.Domain().T1() == 131.892570399495 );

    SECTION( "PointAt(t=0.0)" ) {
        auto point = geometry.PointAt(0.0);

        CHECK( point[0] == Approx(  0) );
        CHECK( point[1] == Approx(-25) );
        CHECK( point[2] == Approx(- 5) );
    }

    SECTION( "PointAt(t=65.9462851997)" ) {
        auto point = geometry.PointAt(65.9462851997);

        REQUIRE( point[0] == Approx( 18.8888888889) );
        REQUIRE( point[1] == Approx(- 6.1111111111) );
        REQUIRE( point[2] == Approx(  4.1666666667) );
    }

    SECTION( "PointAt(t=131.892570399495)" ) {
        auto point = geometry.PointAt(131.892570399495);

        REQUIRE( point[0] == Approx(-25) );
        REQUIRE( point[1] == Approx( 15) );
        REQUIRE( point[2] == Approx(  4) );
    }

    SECTION( "DerivativesAt(t=0.0)" ) {
        auto derivatives = geometry.DerivativesAt(0.0, 3);

        REQUIRE( derivatives[0][0] == Approx(  0) );
        REQUIRE( derivatives[0][1] == Approx(-25) );
        REQUIRE( derivatives[0][2] == Approx(- 5) );

        REQUIRE( derivatives[1][0] == Approx(- 1.81966277  ) );
        REQUIRE( derivatives[1][1] == Approx(  1.2131085134) );
        REQUIRE( derivatives[1][2] == Approx(  0.6065542567) );

        REQUIRE( derivatives[2][0] == Approx(  0.2759310497) );
        REQUIRE( derivatives[2][1] == Approx(- 0.0551862099) );
        REQUIRE( derivatives[2][2] == Approx(- 0.0717420729) );

        REQUIRE( derivatives[3][0] == Approx(- 0.0189682773) );
        REQUIRE( derivatives[3][1] == Approx(  0.0005578905) );
        REQUIRE( derivatives[3][2] == Approx(  0.005523116 ) );
    }

    SECTION( "DerivativesAt(t=65.9462851997)" ) {
        auto derivatives = geometry.DerivativesAt(65.9462851997, 3);

        REQUIRE( derivatives[0][0] == Approx( 18.8888888889) );
        REQUIRE( derivatives[0][1] == Approx(- 6.1111111111) );
        REQUIRE( derivatives[0][2] == Approx(  4.1666666667) );

        REQUIRE( derivatives[1][0] == Approx(  0.2021847522) );
        REQUIRE( derivatives[1][1] == Approx(  0.336974587 ) );
        REQUIRE( derivatives[1][2] == Approx(  0.1010923761) );

        REQUIRE( derivatives[2][0] == Approx(- 0.0122636022) );
        REQUIRE( derivatives[2][1] == Approx(  0.0153295028) );
        REQUIRE( derivatives[2][2] == Approx(- 0.0036790807) );

        REQUIRE( derivatives[3][0] == Approx(- 0.0005578905) );
        REQUIRE( derivatives[3][1] == Approx(- 0.0006508723) );
        REQUIRE( derivatives[3][2] == Approx(  0.0000557891) );
    }

    SECTION( "DerivativesAt(t=131.892570399495)" ) {
        auto derivatives = geometry.DerivativesAt(131.892570399495, 3);

        REQUIRE( derivatives[0][0] == Approx(-25           ) );
        REQUIRE( derivatives[0][1] == Approx( 15           ) );
        REQUIRE( derivatives[0][2] == Approx(  4           ) );

        REQUIRE( derivatives[1][0] == Approx(- 2.4262170267) );
        REQUIRE( derivatives[1][1] == Approx(  2.4262170267) );
        REQUIRE( derivatives[1][2] == Approx(  0.8491759593) );

        REQUIRE( derivatives[2][0] == Approx(- 0.1103724199) );
        REQUIRE( derivatives[2][1] == Approx(  0.3311172597) );
        REQUIRE( derivatives[2][2] == Approx(  0.1269282829) );

        REQUIRE( derivatives[3][0] == Approx(- 0.0044631241) );
        REQUIRE( derivatives[3][1] == Approx(  0.0251050729) );
        REQUIRE( derivatives[3][2] == Approx(  0.0092051934) );
    }
}

TEST_CASE( "Refinement of a BSpline curve", "[CurveGeometry][Refinement][BSpline]" ) {
    CurveGeometry1D curve(2, 4, false);

    curve.SetKnot(0, 0.0);
    curve.SetKnot(1, 0.0);
    curve.SetKnot(2, 1.0);
    curve.SetKnot(3, 2.0);
    curve.SetKnot(4, 2.0);

    curve.SetPole(0, {0.00});
    curve.SetPole(1, {1.25});
    curve.SetPole(2, {3.75});
    curve.SetPole(3, {5.00});

    std::vector<double> newKnots = {0.4, 0.8, 1.2};

    auto result = *curve.Refined(newKnots);

    SECTION( "Check knots" ) {
        CHECK( result.NbKnots() == 8 );

        CHECK( result.Knot(0) == Approx( 0.0 ) );
        CHECK( result.Knot(1) == Approx( 0.0 ) );
        CHECK( result.Knot(2) == Approx( 0.4 ) );
        CHECK( result.Knot(3) == Approx( 0.8 ) );
        CHECK( result.Knot(4) == Approx( 1.0 ) );
        CHECK( result.Knot(5) == Approx( 1.2 ) );
        CHECK( result.Knot(6) == Approx( 2.0 ) );
        CHECK( result.Knot(7) == Approx( 2.0 ) );
    }

    SECTION( "Check poles" ) {
        CHECK( result.NbPoles() == 7 );

        CHECK( result.Pole(0).X() == Approx( 0.00 ) );
        CHECK( result.Pole(1).X() == Approx( 0.50 ) );
        CHECK( result.Pole(2).X() == Approx( 1.50 ) );
        CHECK( result.Pole(3).X() == Approx( 2.25 ) );
        CHECK( result.Pole(4).X() == Approx( 2.75 ) );
        CHECK( result.Pole(5).X() == Approx( 4.00 ) );
        CHECK( result.Pole(6).X() == Approx( 5.00 ) );
    }
}

TEST_CASE( "Refinement of a Nurbs curve", "[CurveGeometry][Refinement][Nurbs]" ) {
    CurveGeometry1D curve(2, 4, true);

    curve.SetKnot(0, 0.0);
    curve.SetKnot(1, 0.0);
    curve.SetKnot(2, 1.0);
    curve.SetKnot(3, 2.0);
    curve.SetKnot(4, 2.0);

    curve.SetPole(0, {0.00});
    curve.SetPole(1, {1.25});
    curve.SetPole(2, {3.75});
    curve.SetPole(3, {5.00});

    curve.SetWeight(0, 1.0);
    curve.SetWeight(1, 5.0);
    curve.SetWeight(2, 1.0);
    curve.SetWeight(3, 1.0);

    std::vector<double> newKnots = {0.4, 0.8, 1.2};

    auto result = *curve.Refined(newKnots);

    SECTION( "Check knots" ) {
        CHECK( result.NbKnots() == 8 );

        CHECK( result.Knot(0) == Approx( 0.0 ) );
        CHECK( result.Knot(1) == Approx( 0.0 ) );
        CHECK( result.Knot(2) == Approx( 0.4 ) );
        CHECK( result.Knot(3) == Approx( 0.8 ) );
        CHECK( result.Knot(4) == Approx( 1.0 ) );
        CHECK( result.Knot(5) == Approx( 1.2 ) );
        CHECK( result.Knot(6) == Approx( 2.0 ) );
        CHECK( result.Knot(7) == Approx( 2.0 ) );
    }

    SECTION( "Check weights" ) {
        CHECK( result.NbPoles() == 7 );

        CHECK( result.Weight(0) == Approx( 1.00 ) );
        CHECK( result.Weight(1) == Approx( 2.60 ) );
        CHECK( result.Weight(2) == Approx( 3.88 ) );
        CHECK( result.Weight(3) == Approx( 3.40 ) );
        CHECK( result.Weight(4) == Approx( 2.60 ) );
        CHECK( result.Weight(5) == Approx( 1.00 ) );
        CHECK( result.Weight(6) == Approx( 1.00 ) );
    }

    SECTION( "Check poles" ) {
        CHECK( result.NbPoles() == 7 );

        CHECK( result.Pole(0).X() == Approx( 0.0000000000 ) );
        CHECK( result.Pole(1).X() == Approx( 0.9615384615 ) );
        CHECK( result.Pole(2).X() == Approx( 1.3144329897 ) );
        CHECK( result.Pole(3).X() == Approx( 1.5441176471 ) );
        CHECK( result.Pole(4).X() == Approx( 1.8269230769 ) );
        CHECK( result.Pole(5).X() == Approx( 4.0000000000 ) );
        CHECK( result.Pole(6).X() == Approx( 5.0000000000 ) );
    }
}
