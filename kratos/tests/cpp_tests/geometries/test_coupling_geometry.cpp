//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License          BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/coupling_geometry.h"

namespace Kratos {
    namespace Testing {

        /// Factory functions
        /** Generates a sample triangle2D3 with Point.
         * Generates a point type right triangle with origin in the origin and leg size 1.
         * @return  Pointer to a triangle2D3
         */
        Triangle2D3<Point>::Pointer GeneratePointTriangle() {
            return Triangle2D3<Point>::Pointer(new Triangle2D3<Point>(
                Point::Pointer(new Point(0.0, 0.0, 0.0)),
                Point::Pointer(new Point(1.0, 0.0, 0.0)),
                Point::Pointer(new Point(0.0, 1.0, 0.0))
                ));
        }

        /** Generates a sample triangle2D3 with Node.
         * @return  Pointer to a triangle2D3
         */
        Triangle2D3<Node<3>>::Pointer GenerateNodeTriangle() {
            return Triangle2D3<Node<3>>::Pointer(new Triangle2D3<Node<3>>(
                Node<3>::Pointer(new Node<3>(1, 1.0, 1.0, 0.0)),
                Node<3>::Pointer(new Node<3>(2, 3.0, 0.5, 0.0)),
                Node<3>::Pointer(new Node<3>(3, 2.5, 2.0, 0.0))
                ));
        }

        /// Tests
        KRATOS_TEST_CASE_IN_SUITE(CouplingNode, KratosCoreGeometriesFastSuite) {
            auto triangle_master = GenerateNodeTriangle();
            auto triangle_slave = GenerateNodeTriangle();
            auto triangle_second_slave = GenerateNodeTriangle();
            auto triangle_third_slave = GenerateNodeTriangle();

            auto coupling_geometry = CouplingGeometry<Node<3>>::Pointer(
                new CouplingGeometry<Node<3>>(triangle_master, triangle_slave));

            KRATOS_CHECK_EQUAL(coupling_geometry->Dimension(), 2);
            KRATOS_CHECK_EQUAL(coupling_geometry->WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(coupling_geometry->LocalSpaceDimension(), 2);

            KRATOS_CHECK_NEAR(coupling_geometry->DomainSize(), 1.375, TOLERANCE);
            KRATOS_CHECK_NEAR(coupling_geometry->Center()[0], 2.1666666667, TOLERANCE);

            std::size_t index = coupling_geometry->AddGeometryPart(triangle_second_slave);
            coupling_geometry->SetGeometryPart(index, triangle_third_slave);

            KRATOS_CHECK_EQUAL(coupling_geometry->NumberOfGeometryParts(), 3);
        }

        KRATOS_TEST_CASE_IN_SUITE(CouplingPoint, KratosCoreGeometriesFastSuite) {
            auto triangle_master = GeneratePointTriangle();
            auto triangle_slave = GeneratePointTriangle();

            auto coupling_geometry = CouplingGeometry<Point>::Pointer(
                new CouplingGeometry<Point>(triangle_master, triangle_slave));

            KRATOS_CHECK_EQUAL(coupling_geometry->NumberOfGeometryParts(), 2);
        }
    } // namespace Testing.
} // namespace Kratos.
