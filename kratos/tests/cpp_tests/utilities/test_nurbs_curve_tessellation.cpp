//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Apostolatos
//					 Tobias Teschemacher
//					 Thomas Oberbichler
//
//
//  Ported from the ANurbs library (https://github.com/oberbichler/ANurbs)
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

#include "containers/pointer_vector.h"

#include "geometries/geometry.h"
#include "geometries/nurbs_curve_geometry.h"

#include "utilities/nurbs_curve_tessellation.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(NurbsCurveTessellation3D,
                          KratosCoreNurbsGeometriesFastSuite)
{
    // Control Points
    PointerVector<Point> points;
    points.push_back(Point::Pointer(new Point(0.0, 0.0)));
    points.push_back(Point::Pointer(new Point(1.0, 1.0)));
    points.push_back(Point::Pointer(new Point(2.0, -1.0)));
    points.push_back(Point::Pointer(new Point(3.0, -1.0)));

    // Knot vector
    Vector knot_vector = ZeroVector(5);
    knot_vector[0] = 0.0;
    knot_vector[1] = 0.0;
    knot_vector[2] = 1.0;
    knot_vector[3] = 2.0;
    knot_vector[4] = 2.0;

    // Polynomial order
    const int p = 2;

    // Create the 2d nurbs curve
    /* const auto geometry = Kratos::make_shared<NurbsCurveGeometry<2, PointerVector<Point>>>(
        points,
        p,
        knot_vector);*/
    const auto geometry = NurbsCurveGeometry<2, PointerVector<Point>>(
        points,
        p,
        knot_vector);

    const auto tessellation = NurbsCurveTessellation<2, PointerVector<Point>>::ComputeTessellation(
        geometry,
        geometry.PolynomialDegree(),
        geometry.DomainInterval(),
        geometry.KnotSpanIntervals(),
        1e-2);

    const double point_expected[17][3] = {{0.0, 0.0, 0.0},
                                    {0.125, 0.2421875, 0.21875},
                                    {0.25, 0.46875, 0.375},
                                    {0.375, 0.6796875, 0.46875},
                                    {0.5, 0.875, 0.5},
                                    {0.625, 1.0546875, 0.46875},
                                    {0.75, 1.21875, 0.375},
                                    {0.875, 1.3671875, 0.21875},
                                    {1.0, 1.5, 0.0},
                                    {1.125, 1.6328125, -0.234375},
                                    {1.25, 1.78125, -0.4375},
                                    {1.375, 1.9453125, -0.609375},
                                    {1.5, 2.125, -0.75},
                                    {1.625, 2.3203125, -0.859375},
                                    {1.75, 2.53125, -0.9375},
                                    {1.875, 2.7578125, -0.984375},
                                    {2.0, 3.0, -1.0}};

    unsigned counter = 0;
    for (auto it = tessellation.begin(); it != tessellation.end(); it++) 
    {
        KRATOS_CHECK_NEAR(it->first, point_expected[counter][0], TOLERANCE);
        KRATOS_CHECK_VECTOR_NEAR(it->second, Point(point_expected[counter][1], point_expected[counter][2]), TOLERANCE)
        counter++;
    }
}

} // namespace Testing.
} // namespace Kratos.
