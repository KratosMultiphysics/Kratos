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
#include "geometries/nurbs_curve_tessellation.h"
#include "nurbs_test_utility.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurveTessellation3D,
        KratosCoreNurbsGeometriesFastSuite)
    {
        NurbsCurveGeometry<2> geometry(
            2,      // degree
            4,      // nb_poles
            false   // IsRational
        );

        geometry.SetKnot(0, 0.0);
        geometry.SetKnot(1, 0.0);
        geometry.SetKnot(2, 1.0);
        geometry.SetKnot(3, 2.0);
        geometry.SetKnot(4, 2.0);
        
        geometry.SetPole(0, NurbsTestUtility::Point2D(0.0,  0.0));
        geometry.SetPole(1, NurbsTestUtility::Point2D(1.0,  1.0));
        geometry.SetPole(2, NurbsTestUtility::Point2D(2.0, -1.0));
        geometry.SetPole(3, NurbsTestUtility::Point2D(3.0, -1.0));

        const auto tessellation =
            NurbsCurveTessellation<NurbsCurveGeometry<2>>::Compute(geometry,
                geometry.GetDomain(), 1e-2);

        NurbsTestUtility::ArrayAlmostEqual(tessellation[0].second,
            NurbsTestUtility::Point2D(0.0000000,  0.000000));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[1].second,
            NurbsTestUtility::Point2D(0.2421875,  0.218750));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[2].second,
            NurbsTestUtility::Point2D(0.4687500,  0.375000));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[3].second,
            NurbsTestUtility::Point2D(0.6796875,  0.468750));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[4].second,
            NurbsTestUtility::Point2D(0.8750000,  0.500000));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[5].second,
            NurbsTestUtility::Point2D(1.0546875,  0.468750));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[6].second,
            NurbsTestUtility::Point2D(1.2187500,  0.375000));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[7].second,
            NurbsTestUtility::Point2D(1.3671875,  0.218750));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[8].second,
            NurbsTestUtility::Point2D(1.5000000,  0.000000));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[9].second,
            NurbsTestUtility::Point2D(1.6328125, -0.234375));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[10].second,
            NurbsTestUtility::Point2D(1.7812500, -0.437500));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[11].second,
            NurbsTestUtility::Point2D(1.9453125, -0.609375));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[12].second,
            NurbsTestUtility::Point2D(2.1250000, -0.750000));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[13].second,
            NurbsTestUtility::Point2D(2.3203125, -0.859375));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[14].second,
            NurbsTestUtility::Point2D(2.5312500, -0.937500));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[15].second,
            NurbsTestUtility::Point2D(2.7578125, -0.984375));
        NurbsTestUtility::ArrayAlmostEqual(tessellation[16].second,
            NurbsTestUtility::Point2D(3.0000000, -1.000000));
    }

} // namespace Testing.
} // namespace Kratos.
