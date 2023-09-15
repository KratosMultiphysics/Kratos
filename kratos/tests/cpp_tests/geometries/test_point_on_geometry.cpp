//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License          BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

#include "geometries/nurbs_curve_geometry.h"
#include "geometries/point_on_geometry.h"

namespace Kratos {
    namespace Testing {

        typedef Node NodeType;

        /// Factory functions
        inline typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer GenerateReferenceCurve3dForPointOnGeometry()
        {
            PointerVector<Node> points(2);

            points(0) = Kratos::make_intrusive<Node>(1, 0, 0, 0);
            points(1) = Kratos::make_intrusive<Node>(2, 10, 0, 0);

            Vector knot_vector = ZeroVector(2);
            knot_vector[0] = 0.0;
            knot_vector[1] = 1.0;

            IndexType p = 1;

            return Kratos::make_shared<NurbsCurveGeometry<3, PointerVector<Node>>>(points, p, knot_vector);
        }

        /// Check dimensions and location
        KRATOS_TEST_CASE_IN_SUITE(PointOnGeometry, KratosCoreGeometriesFastSuite) {
            auto p_nurbs_curve = GenerateReferenceCurve3dForPointOnGeometry();

            array_1d<double, 3> point_coordinates = ZeroVector(3);
            point_coordinates[0] = 0.5;

            PointOnGeometry<PointerVector<Node>, 3, 1> point(point_coordinates, p_nurbs_curve);

            KRATOS_EXPECT_EQ(point.LocalSpaceDimension(), 1);
            KRATOS_EXPECT_EQ(point.WorkingSpaceDimension(), 3);

            /// Check creation of quadrature point geometries
            GeometryType::GeometriesArrayType geometry_vector;
            auto integration_info = point.GetDefaultIntegrationInfo();
            point.CreateQuadraturePointGeometries(geometry_vector, 2, integration_info);

            KRATOS_EXPECT_EQ(geometry_vector[0].LocalSpaceDimension(), 1);
            KRATOS_EXPECT_EQ(geometry_vector[0].WorkingSpaceDimension(), 3);

            std::vector<double> location{ { 5, 0, 0.0 } };
            KRATOS_EXPECT_VECTOR_NEAR(geometry_vector[0].Center(), location, TOLERANCE);
            std::vector<double> shape_functions{ { 0.5, 0.5 } };
            KRATOS_EXPECT_VECTOR_NEAR(row(geometry_vector[0].ShapeFunctionsValues(), 0), shape_functions, TOLERANCE);
        }
    } // namespace Testing.
} // namespace Kratos.
