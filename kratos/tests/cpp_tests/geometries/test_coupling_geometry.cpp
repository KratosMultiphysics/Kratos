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
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/coupling_geometry.h"

namespace Kratos {
    namespace Testing {

        typedef Node<3> NodeType;

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

        Triangle3D3<Point>::Pointer GeneratePointTriangle3D() {
            return Kratos::make_shared<Triangle3D3<Point>>(
                Kratos::make_unique<Point>(0.0, 0.0, 0.0),
                Kratos::make_unique<Point>(1.0, 0.0, 0.0),
                Kratos::make_unique<Point>(0.0, 1.0, 0.0)
                );
        }

        /** Generates a sample triangle2D3 with Node.
         * @return  Pointer to a triangle2D3
         */
        Triangle2D3<NodeType>::Pointer GenerateNodeTriangle() {
            return Kratos::make_shared<Triangle2D3<NodeType>>(
                new NodeType(1, 1.0, 1.0, 0.0),
                new NodeType(2, 3.0, 0.5, 0.0),
                new NodeType(3, 2.5, 2.0, 0.0)
                );
        }

        /// Test with PointType = Node
        KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryNode, KratosCoreGeometriesFastSuite) {
            auto p_triangle_master = GenerateNodeTriangle();
            auto p_triangle_slave = GenerateNodeTriangle();
            auto p_triangle_second_slave = GenerateNodeTriangle();
            auto p_triangle_third_slave = GenerateNodeTriangle();

            auto p_coupling_geometry = Kratos::make_shared<CouplingGeometry<NodeType>>(
                p_triangle_master, p_triangle_slave);

            KRATOS_CHECK_EQUAL(p_coupling_geometry->Dimension(), 2);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->LocalSpaceDimension(), 2);

            KRATOS_CHECK_NEAR(p_coupling_geometry->DomainSize(), 1.375, TOLERANCE);
            KRATOS_CHECK_NEAR(p_coupling_geometry->Center()[0], 2.1666666667, TOLERANCE);

            // Check if master geometry can be found.
            KRATOS_CHECK_NEAR(p_coupling_geometry->GetGeometryPart(CouplingGeometry<Node<3>>::Master).DomainSize(), 1.375, TOLERANCE);
            // Check if geometry 2 can be found.
            KRATOS_CHECK_EQUAL(p_coupling_geometry->GetGeometryPart(CouplingGeometry<Node<3>>::Slave).Dimension(), 2);

            KRATOS_DEBUG_CHECK_EXCEPTION_IS_THROWN(p_coupling_geometry->GetGeometryPart(2), "Index 2 out of range. Composite contains only of: 2 geometries.")

            std::size_t index = p_coupling_geometry->AddGeometryPart(p_triangle_second_slave);
            p_coupling_geometry->SetGeometryPart(index, p_triangle_third_slave);

            KRATOS_CHECK_EQUAL(p_coupling_geometry->NumberOfGeometryParts(), 3);

            // Check the nodes are preserved
            auto& r_geometry_master = *p_triangle_master;
            auto& r_geometry_0 = p_coupling_geometry->GetGeometryPart(0);
            for (std::size_t i = 0; i < r_geometry_master.size(); ++i) {
                KRATOS_CHECK_EQUAL(r_geometry_0[i].Id(), r_geometry_master[i].Id());
                KRATOS_CHECK_LESS(norm_2(r_geometry_0[i].Coordinates() - r_geometry_master[i].Coordinates()), 1.0e-6);
            }
            auto& r_geometry_slave = *p_triangle_slave;
            auto& r_geometry_1 = p_coupling_geometry->GetGeometryPart(1);
            for (std::size_t i = 0; i < r_geometry_slave.size(); ++i) {
                KRATOS_CHECK_EQUAL(r_geometry_1[i].Id(), r_geometry_slave[i].Id());
                KRATOS_CHECK_LESS(norm_2(r_geometry_1[i].Coordinates() - r_geometry_slave[i].Coordinates()), 1.0e-6);
            }
        }

        /// Test with PointType = Point
        KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryPoint, KratosCoreGeometriesFastSuite) {
            auto p_triangle_master = GeneratePointTriangle();
            auto p_triangle_slave = GeneratePointTriangle();

            auto p_coupling_geometry = CouplingGeometry<Point>::Pointer(
                new CouplingGeometry<Point>(p_triangle_master, p_triangle_slave));

            KRATOS_CHECK_EQUAL(p_coupling_geometry->NumberOfGeometryParts(), 2);
        }

        /// Test to exchange the master of the coupling geometry
        KRATOS_TEST_CASE_IN_SUITE(CouplingGeometryChangeMaster, KratosCoreGeometriesFastSuite) {
            auto p_triangle_master = GeneratePointTriangle();
            auto p_triangle_slave = GeneratePointTriangle();

            auto p_triangle_master_2 = GeneratePointTriangle3D();

            auto p_coupling_geometry = CouplingGeometry<Point>::Pointer(
                new CouplingGeometry<Point>(p_triangle_master, p_triangle_slave));

            KRATOS_CHECK_EQUAL(p_coupling_geometry->Dimension(), 2);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->LocalSpaceDimension(), 2);

            p_coupling_geometry->SetGeometryPart(CouplingGeometry<Point>::Master, p_triangle_master_2);

            KRATOS_CHECK_EQUAL(p_coupling_geometry->Dimension(), 2);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->WorkingSpaceDimension(), 3);
            KRATOS_CHECK_EQUAL(p_coupling_geometry->LocalSpaceDimension(), 2);
        }
    } // namespace Testing.
} // namespace Kratos.
