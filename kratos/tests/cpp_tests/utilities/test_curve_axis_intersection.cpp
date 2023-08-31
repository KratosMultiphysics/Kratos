//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License          BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "utilities/curve_axis_intersection.h"

namespace Kratos {
    namespace Testing {

        /// Test Compute Axis Intersection 
        KRATOS_TEST_CASE_IN_SUITE(CurveAxisIntersectionComputeAxisIntersection, KratosCoreFastSuite) {

            /// Generates a sample NurbsCurve:
            PointerVector<Point> points_curve(2);
            points_curve(0) = Kratos::make_shared<Point>(1, -0.23561944901923448, 0);
            points_curve(1) = Kratos::make_shared<Point>(1, -0.11780972450961724, 0);
            Vector knot_vector_curve = ZeroVector(4);
            knot_vector_curve[0] = -0.23561944901923448;
            knot_vector_curve[1] = -0.23561944901923448;
            knot_vector_curve[2] = -0.11780972450961724;
            knot_vector_curve[3] = -0.11780972450961724;
            Vector weights_curve = ZeroVector(2);
            weights_curve[0] = 1;
            weights_curve[1] = 1;
        
            auto p_curve = Kratos::make_shared<
                NurbsCurveGeometry<2, PointerVector<Point>>>(
                    points_curve, 1, knot_vector_curve, weights_curve);

            //Surface span vector
            std::vector<double> surface_spans_u;
            surface_spans_u.resize(3);
            surface_spans_u[0] = 0.0;
            surface_spans_u[1] = 0.5;
            surface_spans_u[2] = 1.0;
            std::vector<double> surface_spans_v;
            surface_spans_v.resize(3);
            surface_spans_v[0] = -0.23561944901923448;
            surface_spans_v[1] = -0.17671458676442586;
            surface_spans_v[2] = -0.11780972450961724;

            std::vector<double> spans;
            CurveAxisIntersection<PointerVector<Point>::value_type>::ComputeAxisIntersection(
                spans,
                *(p_curve.get()), p_curve->DomainInterval().GetT0(), p_curve->DomainInterval().GetT1(),
                surface_spans_u, surface_spans_v,
                1e-6);
                
            KRATOS_EXPECT_EQ(spans.size(), 3);
            KRATOS_EXPECT_NEAR(spans[0], -0.23561944901923448, TOLERANCE);
            KRATOS_EXPECT_NEAR(spans[1], -0.17671458676442586, TOLERANCE);
            KRATOS_EXPECT_NEAR(spans[2], -0.11780972450961724, TOLERANCE);
        }

    } // namespace Testing.
} // namespace Kratos.
