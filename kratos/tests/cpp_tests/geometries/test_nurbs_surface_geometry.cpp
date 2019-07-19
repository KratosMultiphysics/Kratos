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
            2,      // DegreeU
            1,      // DegreeV
            4,      // NbPolesU
            3,      // NbPolesV
            true    // IsRational
        );

        geometry.SetKnotU(0,  0.0);
        geometry.SetKnotU(1,  0.0);
        geometry.SetKnotU(2,  7.5);
        geometry.SetKnotU(3, 15.0);
        geometry.SetKnotU(4, 15.0);

        geometry.SetKnotV(0,  0.0);
        geometry.SetKnotV(1, 10.0);
        geometry.SetKnotV(2, 20.0);

        geometry.SetPole(0, 0, NurbsTestUtility::Point(-10.0, - 5.0, -1.0         ));
        geometry.SetPole(0, 1, NurbsTestUtility::Point(-12.0,   3.0,  3.0         ));
        geometry.SetPole(0, 2, NurbsTestUtility::Point(- 9.0,  11.0, -0.0701928417));
        geometry.SetPole(3,    NurbsTestUtility::Point(- 5.0, - 3.0,  1.0         ));
        geometry.SetPole(1, 1, NurbsTestUtility::Point(- 6.0,   4.0, -2.0         ));
        geometry.SetPole(1, 2, NurbsTestUtility::Point(- 5.0,   7.0,  0.9298071583));
        geometry.SetPole(2, 0, NurbsTestUtility::Point(  0.0, - 4.0, -1.0         ));
        geometry.SetPole(2, 1, NurbsTestUtility::Point(  1.0,   6.0,  5.0         ));
        geometry.SetPole(2, 2, NurbsTestUtility::Point(  0.0,  13.0, -0.2350184214));
        geometry.SetPole(3, 0, NurbsTestUtility::Point(  4.0, - 2.0,  0.0         ));
        geometry.SetPole(3, 1, NurbsTestUtility::Point(  5.0,   4.0, -1.0         ));
        geometry.SetPole(3, 2, NurbsTestUtility::Point(  5.0,  11.0,  0.7649815786));

        geometry.SetWeight(0, 0, 1.0);
        geometry.SetWeight(0, 1, 1.0);
        geometry.SetWeight(0, 2, 1.0);
        geometry.SetWeight(1, 0, 1.0);
        geometry.SetWeight(1, 1, 2.5);
        geometry.SetWeight(1, 2, 1.0);
        geometry.SetWeight(2, 0, 1.0);
        geometry.SetWeight(7,    1.0);
        geometry.SetWeight(2, 2, 1.0);
        geometry.SetWeight(3, 0, 1.0);
        geometry.SetWeight(3, 1, 1.0);
        geometry.SetWeight(3, 2, 1.0);

        KRATOS_CHECK_EQUAL(geometry.GetDimension(), 3);

        KRATOS_CHECK_EQUAL(geometry.GetDegreeU(), 2);
        KRATOS_CHECK_EQUAL(geometry.GetDegreeV(), 1);

        KRATOS_CHECK_EQUAL(geometry.GetNbKnotsU(), 5);
        KRATOS_CHECK_EQUAL(geometry.GetNbKnotsV(), 3);

        KRATOS_CHECK_EQUAL(geometry.GetNbPolesU(), 4);
        KRATOS_CHECK_EQUAL(geometry.GetNbPolesV(), 3);
        KRATOS_CHECK_EQUAL(geometry.GetNbPoles(), 12);

        // GetPointAt

        NurbsTestUtility::ArrayAlmostEqual(geometry.GetPointAt(12, 5),
            NurbsTestUtility::Point(1.0377358491, 1.1320754717, 0.7358490566));

        // GetDerivativesAt

        {
            const auto derivatives = geometry.GetDerivativesAt(12, 5, 2);

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
