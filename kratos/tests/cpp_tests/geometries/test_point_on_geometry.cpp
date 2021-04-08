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

#include "geometries/triangle_3d_3.h"
#include "geometries/point_on_geometry.h"

namespace Kratos {
    namespace Testing {

        typedef Node<3> NodeType;

        /// Factory functions
        inline Triangle3D3<NodeType>::Pointer GenerateNodeTriangleForPointOnGeometry() {
            return Kratos::make_shared<Triangle3D3<NodeType>>(
                new NodeType(4, 0.0, 0.0, 0.0),
                new NodeType(5, 1.0, 0.0, 0.0),
                new NodeType(6, 0, 1.0, 0.0)
                );
        }

        /// Check dimensions and location
        KRATOS_TEST_CASE_IN_SUITE(PointOnGeometry, KratosCoreGeometriesFastSuite) {
            auto p_triangle_master = GenerateNodeTriangleForPointOnGeometry();

            array_1d<double, 3> point_coordinates = ZeroVector(3);
            point_coordinates[0] = 0.5;
            point_coordinates[1] = 0.5;

            PointOnGeometry<PointerVector<Node<3>>, 2, 3> point(point_coordinates, p_triangle_master);

            KRATOS_CHECK_EQUAL(point.Dimension(), 0);
            KRATOS_CHECK_EQUAL(point.LocalSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(point.WorkingSpaceDimension(), 3);

            std::vector<double> location = { 0.33333333333, 0.33333333333, 0.0 };
            KRATOS_CHECK_VECTOR_NEAR(point.Center(), location, TOLERANCE);
        }

    } // namespace Testing.
} // namespace Kratos.
