#include "catch.hpp"

#include <ANurbs/src/Knots.h>

using namespace ANurbs;

TEST_CASE( "Find span", "[Knots]" ) {
    int degree = 3;
    std::vector<double> knots = {0, 0, 0, 1, 2, 2, 3, 3, 3};

    SECTION( "Upper span" ) {
        CHECK( ANurbs::Knots::UpperSpan(degree, knots, 0) == 2 );
        CHECK( ANurbs::Knots::UpperSpan(degree, knots, 1) == 3 );
        CHECK( ANurbs::Knots::UpperSpan(degree, knots, 2) == 5 );
        CHECK( ANurbs::Knots::UpperSpan(degree, knots, 3) == 5 );
    }

    SECTION( "Lower span" ) {
        CHECK( ANurbs::Knots::LowerSpan(degree, knots, 0) == 2 );
        CHECK( ANurbs::Knots::LowerSpan(degree, knots, 1) == 2 );
        CHECK( ANurbs::Knots::LowerSpan(degree, knots, 2) == 3 );
        CHECK( ANurbs::Knots::LowerSpan(degree, knots, 3) == 5 );
    }
}

TEST_CASE( "Knot utilities", "[Knots]" ) {
    int degree = 3;
    int nbKnots = 8;
    int nbPoles = 6;

    SECTION( "Degree" ) {
        CHECK( ANurbs::Knots::Degree(nbKnots, nbPoles) == degree );
    }

    SECTION( "Number of knots" ) {
        CHECK( ANurbs::Knots::NbKnots(degree, nbPoles) == nbKnots );
    }

    SECTION( "Number of poles" ) {
        CHECK( ANurbs::Knots::NbPoles(degree, nbKnots) == nbPoles );
    }

    SECTION( "Number of spans" ) {
        CHECK( ANurbs::Knots::NbSpans(degree, nbKnots) == 3 );
    }
}
