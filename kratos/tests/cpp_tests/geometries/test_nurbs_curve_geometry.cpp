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
#include "geometries/nurbs_curve_geometry.h"
#include "nurbs_test_utility.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveGeometry3D,
        KratosCoreNurbsGeometriesFastSuite)
    {
        NurbsCurveGeometry<3> geometry(
            4,      // degree
            8,      // nb_poles
            true    // is_rational
        );

        geometry.set_knot( 0,   0              );
        geometry.set_knot( 1,   0              );
        geometry.set_knot( 2,   0              );
        geometry.set_knot( 3,   0              );
        geometry.set_knot( 4,  32.9731425998736);
        geometry.set_knot( 5,  65.9462851997473);
        geometry.set_knot( 6,  98.9194277996209);
        geometry.set_knot( 7, 131.892570399495 );
        geometry.set_knot( 8, 131.892570399495 );
        geometry.set_knot( 9, 131.892570399495 );
        geometry.set_knot(10, 131.892570399495 );

        geometry.set_pole(0, NurbsTestUtility::Point(  0, -25, - 5));
        geometry.set_pole(1, NurbsTestUtility::Point(-15, -15,   0));
        geometry.set_pole(2, NurbsTestUtility::Point(  5, - 5, - 3));
        geometry.set_pole(3, NurbsTestUtility::Point( 15, -15,   3));
        geometry.set_pole(4, NurbsTestUtility::Point( 25,   0,   6));
        geometry.set_pole(5, NurbsTestUtility::Point( 15,  15,   6));
        geometry.set_pole(6, NurbsTestUtility::Point(- 5, - 5, - 3));
        geometry.set_pole(7, NurbsTestUtility::Point(-25,  15,   4));

        geometry.set_weight(0, 1.0);
        geometry.set_weight(1, 1.0);
        geometry.set_weight(2, 1.0);
        geometry.set_weight(3, 2.5);
        geometry.set_weight(4, 1.0);
        geometry.set_weight(5, 1.0);
        geometry.set_weight(6, 1.0);
        geometry.set_weight(7, 1.0);

        KRATOS_CHECK_EQUAL(geometry.dimension(), 3);
        KRATOS_CHECK_EQUAL(geometry.is_rational(), true);

        KRATOS_CHECK_EQUAL(geometry.degree(), 4);
        KRATOS_CHECK_EQUAL(geometry.nb_knots(), 11);
        KRATOS_CHECK_EQUAL(geometry.nb_poles(), 8);

        KRATOS_CHECK_EQUAL(geometry.domain().t0(), 0);
        KRATOS_CHECK_EQUAL(geometry.domain().t1(), 131.892570399495);

        // point_at

        NurbsTestUtility::ArrayAlmostEqual(geometry.point_at(0.0),
            NurbsTestUtility::Point(0, -25, -5));
        
        NurbsTestUtility::ArrayAlmostEqual(geometry.point_at(65.9462851997),
            NurbsTestUtility::Point(21.333333, -3.666667, 4.9));
        
        NurbsTestUtility::ArrayAlmostEqual(geometry.point_at(131.892570399495),
            NurbsTestUtility::Point(-25, 15, 4));

        // derivarives_at

        {
            const auto derivatives = geometry.derivatives_at(0.0, 3);

            NurbsTestUtility::ArrayAlmostEqual(derivatives[0],
                NurbsTestUtility::Point(0, -25, -5));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[1],
                NurbsTestUtility::Point(-1.81966277, 1.2131085134, 0.6065542567));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[2],
                NurbsTestUtility::Point(0.2759310497, -0.0551862099, -0.0717420729));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[3],
                NurbsTestUtility::Point(-0.0189682773, 0.0005578905, 0.005523116));
        }

        {
            const auto derivatives = geometry.derivatives_at(65.9462851997, 3);

            NurbsTestUtility::ArrayAlmostEqual(derivatives[0],
                NurbsTestUtility::Point(21.33333333, -3.66666667, 4.9));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[1],
                NurbsTestUtility::Point(0.20218475, 0.33697459, 0.10109238));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[2],
                NurbsTestUtility::Point(-0.0122636, 0.0153295, -0.00367908));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[3],
                NurbsTestUtility::Point(-5.57890509e-04, -6.50872261e-04, 5.57890509e-05));
        }

        {
            const auto derivatives = geometry.derivatives_at(131.892570399495, 3);

            NurbsTestUtility::ArrayAlmostEqual(derivatives[0],
                NurbsTestUtility::Point(-25, 15, 4));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[1],
                NurbsTestUtility::Point(-2.4262170267, 2.4262170267, 0.8491759593));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[2],
                NurbsTestUtility::Point(-0.1103724199, 0.3311172597, 0.1269282829));
            NurbsTestUtility::ArrayAlmostEqual(derivatives[3],
                NurbsTestUtility::Point(-0.0044631241, 0.0251050729, 0.0092051934));
        }
    }

} // namespace Testing.
} // namespace Kratos.
