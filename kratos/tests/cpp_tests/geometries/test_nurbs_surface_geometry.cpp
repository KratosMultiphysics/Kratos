//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/nurbs_surface_geometry.h"
#include "nurbs_test_utility.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfaceGeometry3D,
        KratosCoreNurbsGeometriesFastSuite)
    {
        NurbsSurfaceGeometry<3> geometry(
            2,      // degree u
            1,      // degree v
            4,      // nb_poles u
            3,      // nb_poles v
            true    // is_rational
        );

        geometry.set_knot_u(0,  0.0);
        geometry.set_knot_u(1,  0.0);
        geometry.set_knot_u(2,  7.5);
        geometry.set_knot_u(3, 15.0);
        geometry.set_knot_u(4, 15.0);

        geometry.set_knot_v(0,  0.0);
        geometry.set_knot_v(1, 10.0);
        geometry.set_knot_v(2, 20.0);

        geometry.set_pole(0, 0, NurbsTestUtility::Point(-10.0, - 5.0, -1.0         ));
        geometry.set_pole(0, 1, NurbsTestUtility::Point(-12.0,   3.0,  3.0         ));
        geometry.set_pole(0, 2, NurbsTestUtility::Point(- 9.0,  11.0, -0.0701928417));
        geometry.set_pole(3,    NurbsTestUtility::Point(- 5.0, - 3.0,  1.0         ));
        geometry.set_pole(1, 1, NurbsTestUtility::Point(- 6.0,   4.0, -2.0         ));
        geometry.set_pole(1, 2, NurbsTestUtility::Point(- 5.0,   7.0,  0.9298071583));
        geometry.set_pole(2, 0, NurbsTestUtility::Point(  0.0, - 4.0, -1.0         ));
        geometry.set_pole(2, 1, NurbsTestUtility::Point(  1.0,   6.0,  5.0         ));
        geometry.set_pole(2, 2, NurbsTestUtility::Point(  0.0,  13.0, -0.2350184214));
        geometry.set_pole(3, 0, NurbsTestUtility::Point(  4.0, - 2.0,  0.0         ));
        geometry.set_pole(3, 1, NurbsTestUtility::Point(  5.0,   4.0, -1.0         ));
        geometry.set_pole(3, 2, NurbsTestUtility::Point(  5.0,  11.0,  0.7649815786));

        geometry.set_weight(0, 0, 1.0);
        geometry.set_weight(0, 1, 1.0);
        geometry.set_weight(0, 2, 1.0);
        geometry.set_weight(1, 0, 1.0);
        geometry.set_weight(1, 1, 2.5);
        geometry.set_weight(1, 2, 1.0);
        geometry.set_weight(2, 0, 1.0);
        geometry.set_weight(7,    1.0);
        geometry.set_weight(2, 2, 1.0);
        geometry.set_weight(3, 0, 1.0);
        geometry.set_weight(3, 1, 1.0);
        geometry.set_weight(3, 2, 1.0);

        KRATOS_CHECK_EQUAL(geometry.dimension(), 3);

        KRATOS_CHECK_EQUAL(geometry.degree_u(), 2);
        KRATOS_CHECK_EQUAL(geometry.degree_v(), 1);

        KRATOS_CHECK_EQUAL(geometry.nb_knots_u(), 5);
        KRATOS_CHECK_EQUAL(geometry.nb_knots_v(), 3);

        KRATOS_CHECK_EQUAL(geometry.nb_poles_u(), 4);
        KRATOS_CHECK_EQUAL(geometry.nb_poles_v(), 3);
        KRATOS_CHECK_EQUAL(geometry.nb_poles(), 12);

        // point_at

        NurbsTestUtility::ArrayAlmostEqual(geometry.point_at(12, 5),
            NurbsTestUtility::Point(1.0377358491, 1.1320754717, 0.7358490566));

        // derivarives_at

        {
            const auto derivatives = geometry.derivatives_at(12, 5, 2);

            NurbsTestUtility::ArrayAlmostEqual(derivatives[0],
                NurbsTestUtility::Point( 1.0377358491,  1.1320754717,  0.7358490566));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[1],
                NurbsTestUtility::Point( 1.1712353151, -0.0830663344, -0.1483327400));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[2],
                NurbsTestUtility::Point(-0.0004271983,  0.8173727305,  0.2294054824));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[3],
                NurbsTestUtility::Point( 0.0334131158,  0.0214196208, -0.1713942830));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[4],
                NurbsTestUtility::Point( 0.0499024922, -0.0351431047, -0.0293933919));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[5],
                NurbsTestUtility::Point( 9.6724e-6   , -0.0185065524, -0.0051940864));
        }
    }

} // namespace Testing.
} // namespace Kratos.
