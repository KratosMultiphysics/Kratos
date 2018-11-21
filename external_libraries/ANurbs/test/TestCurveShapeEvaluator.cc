#include "catch.hpp"

#include <ANurbs/src/CurveShapeEvaluator.h>

#include <functional>

using namespace ANurbs;


TEST_CASE( "Shape functions of a curve", "[CurveShapeEvaluator]" ) {
    int degree = 2;
    int order = 2;

    CurveShapeEvaluator<double> shape(degree, order);

    SECTION( "Compute values of a B-Spline" ) {
        std::vector<double> knots = {1, 1, 3, 3};

        double t = 2;

        shape.Compute(knots, t);

        CHECK( shape(0, 0) == Approx( 0.25) );
        CHECK( shape(0, 1) == Approx( 0.50) );
        CHECK( shape(0, 2) == Approx( 0.25) );
        
        CHECK( shape(1, 0) == Approx(-0.50) );
        CHECK( shape(1, 1) == Approx( 0.00) );
        CHECK( shape(1, 2) == Approx( 0.50) );
        
        CHECK( shape(2, 0) == Approx( 0.50) );
        CHECK( shape(2, 1) == Approx(-1.00) );
        CHECK( shape(2, 2) == Approx( 0.50) );
    }

    SECTION( "Compute values of a equal weighted Nurbs" ) {
        std::vector<double> knots = {1, 1, 3, 3};
        std::vector<double> weights = {2, 2, 2};

        double t = 2;

        shape.Compute(knots, weights, t);

        CHECK( shape(0, 0) == Approx( 0.25) );
        CHECK( shape(0, 1) == Approx( 0.50) );
        CHECK( shape(0, 2) == Approx( 0.25) );
        
        CHECK( shape(1, 0) == Approx(-0.50) );
        CHECK( shape(1, 1) == Approx( 0.00) );
        CHECK( shape(1, 2) == Approx( 0.50) );
        
        CHECK( shape(2, 0) == Approx( 0.50) );
        CHECK( shape(2, 1) == Approx(-1.00) );
        CHECK( shape(2, 2) == Approx( 0.50) );
    }

    SECTION( "Compute values of a weighted Nurbs" ) {
        std::vector<double> knots = {1, 1, 3, 3};
        std::vector<double> weights = {1, 1.5, 1};

        double t = 2;

        shape.Compute(knots, weights, t);

        CHECK( shape(0, 0) == Approx( 0.20) );
        CHECK( shape(0, 1) == Approx( 0.60) );
        CHECK( shape(0, 2) == Approx( 0.20) );
        
        CHECK( shape(1, 0) == Approx(-0.50) );
        CHECK( shape(1, 1) == Approx( 0.00) );
        CHECK( shape(1, 2) == Approx( 0.50) );
        
        CHECK( shape(2, 0) == Approx( 0.50) );
        CHECK( shape(2, 1) == Approx(-1.00) );
        CHECK( shape(2, 2) == Approx( 0.50) );
    }
}

TEST_CASE( "Shape functions of a curve 2", "[CurveShapeEvaluator]" ) {
    int degree = 1;
    int order = 1;

    CurveShapeEvaluator<double> shape(degree, order);

    std::vector<double> knots = {0, 5};

    double t;
    
    t = 1.0566243270259357;

    shape.Compute(knots, t);

    CHECK( shape(0, 0) == Approx( 0.788675134594813  ) );
    CHECK( shape(0, 1) == Approx( 0.21132486540518713) );
    
    CHECK( shape(1, 0) == Approx(-0.2                ) );
    CHECK( shape(1, 1) == Approx( 0.2                ) );

    t = 3.9433756729740645;

    shape.Compute(knots, t);

    CHECK( shape(0, 0) == Approx( 0.2113248654051871 ) );
    CHECK( shape(0, 1) == Approx( 0.788675134594813  ) );
    
    CHECK( shape(1, 0) == Approx(-0.2                ) );
    CHECK( shape(1, 1) == Approx( 0.2                ) );
}
