#include "catch.hpp"

#include <ANurbs/src/Curve.h>
#include <ANurbs/src/CurveGeometry.h>
#include <ANurbs/src/CurveTessellation.h>

using namespace ANurbs;


TEST_CASE( "Tessellate a planar curve", "[CurveTessellation][Curve][2D]" ) {
    // create curve geometry

    int degree = 2;
    int nbPoles = 4;
    bool isRational = false;

    Pointer<CurveGeometry2D> geometry = New<CurveGeometry2D>(degree, nbPoles, isRational);

    geometry->SetKnot(0, 0.0);
    geometry->SetKnot(1, 0.0);
    geometry->SetKnot(2, 1.0);
    geometry->SetKnot(3, 2.0);
    geometry->SetKnot(4, 2.0);

    geometry->SetPole(0, {0.0,  0.0});
    geometry->SetPole(1, {1.0,  1.0});
    geometry->SetPole(2, {2.0, -1.0});
    geometry->SetPole(3, {3.0, -1.0});

    // create curve

    Curve2D curve(geometry);

    // tessellate

    const double flatness = 1e-2;

    CurveTessellation2D tessellation;

    tessellation.Compute(curve, flatness);


    SECTION( "Check points" ) {
        // FIXME: look for a better way to test the result

        CHECK( tessellation.NbPoints() == 17 );

        CHECK( tessellation.Point( 0).Y() == Approx( 0.0      ) );
        CHECK( tessellation.Point( 0).X() == Approx( 0.0      ) );

        CHECK( tessellation.Point( 1).Y() == Approx( 0.21875  ) );
        CHECK( tessellation.Point( 1).X() == Approx( 0.2421875) );

        CHECK( tessellation.Point( 2).Y() == Approx( 0.375    ) );
        CHECK( tessellation.Point( 2).X() == Approx( 0.46875  ) );

        CHECK( tessellation.Point( 3).Y() == Approx( 0.46875  ) );
        CHECK( tessellation.Point( 3).X() == Approx( 0.6796875) );

        CHECK( tessellation.Point( 4).Y() == Approx( 0.5      ) );
        CHECK( tessellation.Point( 4).X() == Approx( 0.875    ) );

        CHECK( tessellation.Point( 5).Y() == Approx( 0.46875  ) );
        CHECK( tessellation.Point( 5).X() == Approx( 1.0546875) );

        CHECK( tessellation.Point( 6).Y() == Approx( 0.375    ) );
        CHECK( tessellation.Point( 6).X() == Approx( 1.21875  ) );

        CHECK( tessellation.Point( 7).Y() == Approx( 0.21875  ) );
        CHECK( tessellation.Point( 7).X() == Approx( 1.3671875) );

        CHECK( tessellation.Point( 8).Y() == Approx( 0.0      ) );
        CHECK( tessellation.Point( 8).X() == Approx( 1.5      ) );

        CHECK( tessellation.Point( 9).Y() == Approx(-0.234375 ) );
        CHECK( tessellation.Point( 9).X() == Approx( 1.6328125) );

        CHECK( tessellation.Point(10).Y() == Approx(-0.4375   ) );
        CHECK( tessellation.Point(10).X() == Approx( 1.78125  ) );

        CHECK( tessellation.Point(11).Y() == Approx(-0.609375 ) );
        CHECK( tessellation.Point(11).X() == Approx( 1.9453125) );

        CHECK( tessellation.Point(12).Y() == Approx(-0.75     ) );
        CHECK( tessellation.Point(12).X() == Approx( 2.125    ) );

        CHECK( tessellation.Point(13).Y() == Approx(-0.859375 ) );
        CHECK( tessellation.Point(13).X() == Approx( 2.3203125) );

        CHECK( tessellation.Point(14).Y() == Approx(-0.9375   ) );
        CHECK( tessellation.Point(14).X() == Approx( 2.53125  ) );

        CHECK( tessellation.Point(15).Y() == Approx(-0.984375 ) );
        CHECK( tessellation.Point(15).X() == Approx( 2.7578125) );

        CHECK( tessellation.Point(16).Y() == Approx(-1.0      ) );
        CHECK( tessellation.Point(16).X() == Approx( 3.0      ) );
    }
}

TEST_CASE( "Tessellate a spatial curve", "[CurveTessellation][Curve][3D]" ) {
    int degree = 2;
    int nbPoles = 4;
    bool isRational = false;

    Pointer<CurveGeometry3D> geometry = New<CurveGeometry3D>(degree, nbPoles,
        isRational);

    geometry->SetKnot(0, 0.0);
    geometry->SetKnot(1, 0.0);
    geometry->SetKnot(2, 1.0);
    geometry->SetKnot(3, 2.0);
    geometry->SetKnot(4, 2.0);

    geometry->SetPole(0, {0.0,  0.0, -1.0});
    geometry->SetPole(1, {1.0,  1.0,  1.0});
    geometry->SetPole(2, {2.0, -1.0,  0.0});
    geometry->SetPole(3, {3.0, -1.0,  0.0});

    Curve3D curve(geometry);

    const double flatness = 1e-2;

    CurveTessellation3D tessellation;

    tessellation.Compute(curve, flatness);

    SECTION( "Check points" ) {
        // FIXME: look for a better way to test the result

        CHECK( tessellation.NbPoints() == 20 );

        CHECK( tessellation.Point( 0).X() == Approx( 0.0        ) );
        CHECK( tessellation.Point( 0).Y() == Approx( 0.0        ) );
        CHECK( tessellation.Point( 0).Z() == Approx(-1.0        ) );

        CHECK( tessellation.Point( 1).X() == Approx( 0.2421875  ) );
        CHECK( tessellation.Point( 1).Y() == Approx( 0.21875    ) );
        CHECK( tessellation.Point( 1).Z() == Approx(-0.5390625  ) );

        CHECK( tessellation.Point( 2).X() == Approx( 0.46875    ) );
        CHECK( tessellation.Point( 2).Y() == Approx( 0.375      ) );
        CHECK( tessellation.Point( 2).Z() == Approx(-0.15625    ) );

        CHECK( tessellation.Point( 3).X() == Approx( 0.6796875  ) );
        CHECK( tessellation.Point( 3).Y() == Approx( 0.46875    ) );
        CHECK( tessellation.Point( 3).Z() == Approx( 0.1484375  ) );

        CHECK( tessellation.Point( 4).X() == Approx( 0.875      ) );
        CHECK( tessellation.Point( 4).Y() == Approx( 0.5        ) );
        CHECK( tessellation.Point( 4).Z() == Approx( 0.375      ) );

        CHECK( tessellation.Point( 5).X() == Approx( 0.966796875) );
        CHECK( tessellation.Point( 5).Y() == Approx( 0.4921875  ) );
        CHECK( tessellation.Point( 5).Z() == Approx( 0.458984375) );

        CHECK( tessellation.Point( 6).X() == Approx( 1.0546875  ) );
        CHECK( tessellation.Point( 6).Y() == Approx( 0.46875    ) );
        CHECK( tessellation.Point( 6).Z() == Approx( 0.5234375  ) );

        CHECK( tessellation.Point( 7).X() == Approx( 1.138671875) );
        CHECK( tessellation.Point( 7).Y() == Approx( 0.4296875  ) );
        CHECK( tessellation.Point( 7).Z() == Approx( 0.568359375) );

        CHECK( tessellation.Point( 8).X() == Approx( 1.21875    ) );
        CHECK( tessellation.Point( 8).Y() == Approx( 0.375      ) );
        CHECK( tessellation.Point( 8).Z() == Approx( 0.59375    ) );

        CHECK( tessellation.Point( 9).X() == Approx( 1.294921875) );
        CHECK( tessellation.Point( 9).Y() == Approx( 0.3046875  ) );
        CHECK( tessellation.Point( 9).Z() == Approx( 0.599609375) );

        CHECK( tessellation.Point(10).X() == Approx( 1.3671875  ) );
        CHECK( tessellation.Point(10).Y() == Approx( 0.21875    ) );
        CHECK( tessellation.Point(10).Z() == Approx( 0.5859375  ) );

        CHECK( tessellation.Point(11).X() == Approx( 1.5        ) );
        CHECK( tessellation.Point(11).Y() == Approx( 0.0        ) );
        CHECK( tessellation.Point(11).Z() == Approx( 0.5        ) );

        CHECK( tessellation.Point(12).X() == Approx( 1.6328125  ) );
        CHECK( tessellation.Point(12).Y() == Approx(-0.234375   ) );
        CHECK( tessellation.Point(12).Z() == Approx( 0.3828125  ) );

        CHECK( tessellation.Point(13).X() == Approx( 1.78125    ) );
        CHECK( tessellation.Point(13).Y() == Approx(-0.4375     ) );
        CHECK( tessellation.Point(13).Z() == Approx( 0.28125    ) );

        CHECK( tessellation.Point(14).X() == Approx( 1.9453125  ) );
        CHECK( tessellation.Point(14).Y() == Approx(-0.609375   ) );
        CHECK( tessellation.Point(14).Z() == Approx( 0.1953125  ) );

        CHECK( tessellation.Point(15).X() == Approx( 2.125      ) );
        CHECK( tessellation.Point(15).Y() == Approx(-0.75       ) );
        CHECK( tessellation.Point(15).Z() == Approx( 0.125      ) );

        CHECK( tessellation.Point(16).X() == Approx( 2.3203125  ) );
        CHECK( tessellation.Point(16).Y() == Approx(-0.859375   ) );
        CHECK( tessellation.Point(16).Z() == Approx( 0.0703125  ) );

        CHECK( tessellation.Point(17).X() == Approx( 2.53125    ) );
        CHECK( tessellation.Point(17).Y() == Approx(-0.9375     ) );
        CHECK( tessellation.Point(17).Z() == Approx( 0.03125    ) );

        CHECK( tessellation.Point(18).X() == Approx( 2.7578125  ) );
        CHECK( tessellation.Point(18).Y() == Approx(-0.984375   ) );
        CHECK( tessellation.Point(18).Z() == Approx( 0.0078125  ) );

        CHECK( tessellation.Point(19).X() == Approx( 3.0        ) );
        CHECK( tessellation.Point(19).Y() == Approx(-1.0        ) );
        CHECK( tessellation.Point(19).Z() == Approx( 0.0        ) );
    }
}
