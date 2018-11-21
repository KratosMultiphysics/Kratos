#include "catch.hpp"

#include <ANurbs/src/Grid.h>
#include <ANurbs/src/SurfaceShapeEvaluator.h>

using namespace ANurbs;


TEST_CASE( "Shape functions of a surface", "[SurfaceShapeEvaluator]" ) {
    int degreeU = 2;
    int degreeV = 1;
    int order = 2;

    SurfaceShapeEvaluator<double> shape(degreeU, degreeV, order);

    SECTION( "Compute values of a B-Spline" ) {
        std::vector<double> knotsU = {0, 0, 7.5, 15, 15.0};
        std::vector<double> knotsV = {0, 10, 20};

        double u = 12;
        double v = 5;

        shape.Compute(knotsU, knotsV, u, v);

        REQUIRE( shape(0, 0, 0) == Approx( 0.04        ) );
        REQUIRE( shape(0, 0, 1) == Approx( 0.04        ) );
        REQUIRE( shape(0, 1, 0) == Approx( 0.28        ) );
        REQUIRE( shape(0, 1, 1) == Approx( 0.28        ) );
        REQUIRE( shape(0, 2, 0) == Approx( 0.18        ) );
        REQUIRE( shape(0, 2, 1) == Approx( 0.18        ) );

        REQUIRE( shape(1, 0, 0) == Approx(-0.0266666667) );
        REQUIRE( shape(1, 0, 1) == Approx(-0.0266666667) );
        REQUIRE( shape(1, 1, 0) == Approx(-0.0533333333) );
        REQUIRE( shape(1, 1, 1) == Approx(-0.0533333333) );
        REQUIRE( shape(1, 2, 0) == Approx( 0.08        ) );
        REQUIRE( shape(1, 2, 1) == Approx( 0.08        ) );
    }

    SECTION( "Compute values of a equal weighted Nurbs" ) {
        std::vector<double> knotsU = {0, 0, 7.5, 15, 15.0};
        std::vector<double> knotsV = {0, 10, 20};

        Grid<double> weights(4, 3);
        weights.SetValue(0, 0, 2.0);
        weights.SetValue(0, 1, 2.0);
        weights.SetValue(0, 2, 2.0);
        weights.SetValue(1, 0, 2.0);
        weights.SetValue(1, 1, 2.0);
        weights.SetValue(1, 2, 2.0);
        weights.SetValue(2, 0, 2.0);
        weights.SetValue(2, 1, 2.0);
        weights.SetValue(2, 2, 2.0);
        weights.SetValue(3, 0, 2.0);
        weights.SetValue(3, 1, 2.0);
        weights.SetValue(3, 2, 2.0);

        double u = 12;
        double v = 5;

        shape.Compute(knotsU, knotsV, weights, u, v);

        REQUIRE( shape(0, 0, 0) == Approx( 0.04        ) );
        REQUIRE( shape(0, 0, 1) == Approx( 0.04        ) );
        REQUIRE( shape(0, 1, 0) == Approx( 0.28        ) );
        REQUIRE( shape(0, 1, 1) == Approx( 0.28        ) );
        REQUIRE( shape(0, 2, 0) == Approx( 0.18        ) );
        REQUIRE( shape(0, 2, 1) == Approx( 0.18        ) );

        REQUIRE( shape(1, 0, 0) == Approx(-0.0266666667) );
        REQUIRE( shape(1, 0, 1) == Approx(-0.0266666667) );
        REQUIRE( shape(1, 1, 0) == Approx(-0.0533333333) );
        REQUIRE( shape(1, 1, 1) == Approx(-0.0533333333) );
        REQUIRE( shape(1, 2, 0) == Approx( 0.08        ) );
        REQUIRE( shape(1, 2, 1) == Approx( 0.08        ) );
    }

    SECTION( "Compute values of a weighted Nurbs" ) {
        std::vector<double> knotsU = {0, 0, 7.5, 15, 15.0};
        std::vector<double> knotsV = {0, 10, 20};

        Grid<double> weights(4, 3);
        weights.SetValue(0, 0, 1.0);
        weights.SetValue(0, 1, 1.0);
        weights.SetValue(0, 2, 1.0);
        weights.SetValue(1, 0, 1.0);
        weights.SetValue(1, 1, 2.5); // <<<
        weights.SetValue(1, 2, 1.0);
        weights.SetValue(2, 0, 1.0);
        weights.SetValue(2, 1, 1.0);
        weights.SetValue(2, 2, 1.0);
        weights.SetValue(3, 0, 1.0);
        weights.SetValue(3, 1, 1.0);
        weights.SetValue(3, 2, 1.0);

        double u = 12;
        double v = 5;

        shape.Compute(knotsU, knotsV, weights, u, v);

        REQUIRE( shape(0, 0, 0) == Approx( 0.0377358491) );
        REQUIRE( shape(0, 0, 1) == Approx( 0.0943396226) );
        REQUIRE( shape(0, 1, 0) == Approx( 0.2641509434) );
        REQUIRE( shape(0, 1, 1) == Approx( 0.2641509434) );
        REQUIRE( shape(0, 2, 0) == Approx( 0.1698113208) );
        REQUIRE( shape(0, 2, 1) == Approx( 0.1698113208) );

        REQUIRE( shape(1, 0, 0) == Approx(-0.0237332384) );
        REQUIRE( shape(1, 0, 1) == Approx(-0.0593330960) );
        REQUIRE( shape(1, 1, 0) == Approx(-0.0403465053) );
        REQUIRE( shape(1, 1, 1) == Approx(-0.0403465053) );
        REQUIRE( shape(1, 2, 0) == Approx( 0.0818796725) );
        REQUIRE( shape(1, 2, 1) == Approx( 0.0818796725) );
    }
}
