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
            return Kratos::make_shared<Triangle2D3<Point>>(
                Kratos::make_unique<Point>(0.0, 0.0, 0.0),
                Kratos::make_unique<Point>(1.0, 0.0, 0.0),
                Kratos::make_unique<Point>(0.0, 1.0, 0.0)
                );
        }

        /** Generates a sample triangle2D3 with Node.
         * @return  Pointer to a triangle2D3
         */
        Triangle2D3<Node<3>>::Pointer GenerateNodeTriangle() {
            return Kratos::make_shared<Triangle2D3<Node<3>>>(
                new Node<3>(1, 1.0, 1.0, 0.0),
                new Node<3>(2, 3.0, 0.5, 0.0),
                new Node<3>(3, 2.5, 2.0, 0.0)
                );
        }

        /// Tests
        KRATOS_TEST_CASE_IN_SUITE(CouplingNode, KratosCoreGeometriesFastSuite) {
            auto p_triangle_master = GenerateNodeTriangle();
            auto p_triangle_slave = GenerateNodeTriangle();
            auto p_triangle_second_slave = GenerateNodeTriangle();
            auto p_triangle_third_slave = GenerateNodeTriangle();

            auto p_coupling_geometry = Kratos::make_shared<CouplingGeometry<Node<3>>>(
                p_triangle_master, p_triangle_slave);

            KRATOS_CHECK_EQUAL(p_coupling_geometry->Dimension(), 2);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->LocalSpaceDimension(), 2);

            KRATOS_CHECK_NEAR(p_coupling_geometry->DomainSize(), 1.375, TOLERANCE);
            KRATOS_CHECK_NEAR(p_coupling_geometry->Center()[0], 2.1666666667, TOLERANCE);

            //check if master geometry can be found.
            KRATOS_CHECK_NEAR(p_coupling_geometry->GetGeometryPart(CouplingGeometry<Node<3>>::Master).DomainSize(), 1.375, TOLERANCE);
            //check if geometry 2 can be found.
            KRATOS_CHECK_EQUAL(p_coupling_geometry->GetGeometryPart(CouplingGeometry<Node<3>>::Slave).Dimension(), 2);

            KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(p_coupling_geometry->GetGeometryPart(2), "Index 2 out of range. Composite contains only of: 2 geometries.")

            std::size_t index = p_coupling_geometry->AddGeometryPart(p_triangle_second_slave);
            p_coupling_geometry->SetGeometryPart(index, p_triangle_third_slave);

            KRATOS_CHECK_EQUAL(p_coupling_geometry->NumberOfGeometryParts(), 3);
        }

        KRATOS_TEST_CASE_IN_SUITE(CouplingPoint, KratosCoreGeometriesFastSuite) {
            auto p_triangle_master = GeneratePointTriangle();
            auto p_triangle_slave = GeneratePointTriangle();

            auto p_coupling_geometry = CouplingGeometry<Point>::Pointer(
                new CouplingGeometry<Point>(p_triangle_master, p_triangle_slave));

            KRATOS_CHECK_EQUAL(p_coupling_geometry->NumberOfGeometryParts(), 2);
        }
    } // namespace Testing.
} // namespace Kratos.
