#include "catch.hpp"

#include <ANurbs/src/Point.h>

using namespace ANurbs;


TEST_CASE( "Create a Point2D by initializer list", "[Point]" ) {
    Point<double, 2> point = {1.2, 3.4};

    SECTION( "Access values by index" ) {
        CHECK( point[0] == 1.2 );
        CHECK( point[1] == 3.4 );
    }

    SECTION( "Access values by XY" ) {
        CHECK( point.X() == 1.2 );
        CHECK( point.Y() == 3.4 );
    }
}

TEST_CASE( "Assign Point2D by index", "[Point]" ) {
    Point<double, 2> point;
    point[0] = 1.2;
    point[1] = 3.4;

    SECTION( "Access values by index" ) {
        CHECK( point[0] == 1.2 );
        CHECK( point[1] == 3.4 );
    }

    SECTION( "Access values by XY" ) {
        CHECK( point.X() == 1.2 );
        CHECK( point.Y() == 3.4 );
    }
}

TEST_CASE( "Assign Point2D by SetX/Y", "[Point]" ) {
    Point<double, 2> point;
    point.SetX(1.2);
    point.SetY(3.4);

    SECTION( "Access values by index" ) {
        CHECK( point[0] == 1.2 );
        CHECK( point[1] == 3.4 );
    }

    SECTION( "Access values by XY" ) {
        CHECK( point.X() == 1.2 );
        CHECK( point.Y() == 3.4 );
    }
}

TEST_CASE( "Point3D by initializer list", "[Point]" ) {
    Point<double, 3> point = {1.2, 3.4, 5.6};

    SECTION( "Access values by index" ) {
        CHECK( point[0] == 1.2 );
        CHECK( point[1] == 3.4 );
        CHECK( point[2] == 5.6 );
    }

    SECTION( "Access values by XYZ" ) {
        CHECK( point.X() == 1.2 );
        CHECK( point.Y() == 3.4 );
        CHECK( point.Z() == 5.6 );
    }
}

TEST_CASE( "Assign Point3D by index", "[Point]" ) {
    Point<double, 3> point;
    point[0] = 1.2;
    point[1] = 3.4;
    point[2] = 5.6;

    SECTION( "Access values by index" ) {
        CHECK( point[0] == 1.2 );
        CHECK( point[1] == 3.4 );
        CHECK( point[2] == 5.6 );
    }

    SECTION( "Access values by XYZ" ) {
        CHECK( point.X() == 1.2 );
        CHECK( point.Y() == 3.4 );
        CHECK( point.Z() == 5.6 );
    }
}

TEST_CASE( "Assign Point3D by SetX/Y/Z", "[Point]" ) {
    Point<double, 3> point;
    point.SetX(1.2);
    point.SetY(3.4);
    point.SetZ(5.6);

    SECTION( "Access values by index" ) {
        CHECK( point[0] == 1.2 );
        CHECK( point[1] == 3.4 );
        CHECK( point[2] == 5.6 );
    }

    SECTION( "Access values by XYZ" ) {
        CHECK( point.X() == 1.2 );
        CHECK( point.Y() == 3.4 );
        CHECK( point.Z() == 5.6 );
    }
}
