//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:      BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors: Philipp Bucher
//
//

// Project includes
#include "testing/testing.h"
#include "includes/serializer.h"
#include "includes/file_serializer.h"
#include "includes/stream_serializer.h"

#include "utilities/quadrature_points_utility.h"

#include "geometries/point.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/quadrature_point_geometry.h"
#include "geometries/geometry_shape_function_container.h"

namespace Kratos {
    namespace Testing {
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Node<3> NodeType;
        typedef Geometry<Node<3>> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;

        typedef QuadraturePointGeometry<NodeType, 2, 2> QPGeometryType;
        typedef typename QPGeometryType::Pointer QPGeometryPointerType;

        /** Generates a point type sample triangle2D3.
         * Generates a point type right triangle with origin in the origin and leg size 1.
         * @return  Pointer to a triangle2D3
         */
        GeometryPointerType GeneratePointsTriangle2D3() {
            return Kratos::make_shared<Triangle2D3<NodeType>>(
                new NodeType(1, 0.0, 0.0, 0.0),
                new NodeType(2, 1.0, 0.0, 0.0),
                new NodeType(3, 0.0, 1.0, 0.0)
                );
        }

        GeometryPointerType GenerateQuadraturePointGeometry() {
            auto triangle = GeneratePointsTriangle2D3();

            //auto integration_method = triangle->GetDefaultIntegrationMethod();
            auto integration_points = triangle->IntegrationPoints();

            auto r_N = triangle->ShapeFunctionsValues();

            Matrix N_i = ZeroMatrix(1, triangle->size());
            for (IndexType j = 0; j < triangle->size(); ++j) {
                N_i(0, j) = r_N(0, j);
            }
            Matrix DN_De = triangle->ShapeFunctionLocalGradient(0);

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                GeometryData::GI_GAUSS_1,
                integration_points[0],
                N_i,
                DN_De);

            GeometryPointerType p_this_quadrature_point(
                Kratos::make_shared<QuadraturePointGeometry<Node<3>, 2, 2>>(
                    triangle->Points(),
                    data_container,
                    triangle.get()));

            return p_this_quadrature_point;
        }

        QPGeometryPointerType GenerateQuadraturePointGeometry2() {
            auto triangle = GeneratePointsTriangle2D3();

            auto integration_points = triangle->IntegrationPoints();

            auto r_N = triangle->ShapeFunctionsValues();

            Matrix N_i = ZeroMatrix(1, triangle->size());
            for (IndexType j = 0; j < triangle->size(); ++j) {
                N_i(0, j) = r_N(0, j);
            }
            const Matrix& DN_De = triangle->ShapeFunctionLocalGradient(0);
            DenseVector<Matrix> derivatives(1);
            derivatives[0] = DN_De;

            return Kratos::make_shared<QuadraturePointGeometry<Node<3>, 2, 2>>(
                triangle->Points(),
                integration_points[0],
                N_i,
                derivatives,
                triangle.get());
        }

        KRATOS_TEST_CASE_IN_SUITE(QuadraturePointUtility, KratosCoreFastSuite)
        {
            auto triangle = GeneratePointsTriangle2D3();

            auto quadrature_points = CreateQuadraturePointsUtility<NodeType>::Create(
                triangle, GeometryData::GI_GAUSS_3);

            KRATOS_CHECK_EQUAL(quadrature_points.size(), 4);
            KRATOS_CHECK_EQUAL(quadrature_points[0]->size(), 3);
            KRATOS_CHECK_EQUAL(quadrature_points[0]->WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(quadrature_points[0]->LocalSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(quadrature_points[0]->Dimension(), 2);

            KRATOS_CHECK_MATRIX_NEAR(quadrature_points[0]->ShapeFunctionsLocalGradients()[0], triangle->ShapeFunctionsLocalGradients(GeometryData::GI_GAUSS_3)[0], 1e-6);
            KRATOS_CHECK_MATRIX_NEAR(quadrature_points[1]->ShapeFunctionsLocalGradients()[0], triangle->ShapeFunctionsLocalGradients(GeometryData::GI_GAUSS_3)[1], 1e-6);
        }

        KRATOS_TEST_CASE_IN_SUITE(QuadraturePointGeometry2d, KratosCoreFastSuite)
        {
            auto p_this_quadrature_point = GenerateQuadraturePointGeometry();

            KRATOS_CHECK_EQUAL(p_this_quadrature_point->size(), 3);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->LocalSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->Dimension(), 2);

            KRATOS_CHECK_EQUAL(p_this_quadrature_point->ShapeFunctionsValues().size1(), 1);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->ShapeFunctionsValues().size2(), 3);
        }

        KRATOS_TEST_CASE_IN_SUITE(QuadraturePointGeometry2dCopyConstructor, KratosCoreFastSuite)
        {
            QPGeometryPointerType p_this_quadrature_point = GenerateQuadraturePointGeometry2();

            QPGeometryType geom(*p_this_quadrature_point);

            KRATOS_CHECK_EQUAL(geom.size(), 3);
            KRATOS_CHECK_EQUAL(geom.WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(geom.LocalSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(geom.Dimension(), 2);

            KRATOS_CHECK_EQUAL(geom.ShapeFunctionsValues().size1(), 1);
            KRATOS_CHECK_EQUAL(geom.ShapeFunctionsValues().size2(), 3);

            KRATOS_CHECK_MATRIX_NEAR(
                geom.ShapeFunctionsValues(),
                p_this_quadrature_point->ShapeFunctionsValues(),
                1e-6);
        }

        KRATOS_TEST_CASE_IN_SUITE(QuadraturePointGeometry2dAssignmentOperator, KratosCoreFastSuite)
        {
            QPGeometryPointerType p_this_quadrature_point = GenerateQuadraturePointGeometry2();

            QPGeometryType geom = *p_this_quadrature_point;

            KRATOS_CHECK_EQUAL(geom.size(), 3);
            KRATOS_CHECK_EQUAL(geom.WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(geom.LocalSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(geom.Dimension(), 2);

            KRATOS_CHECK_EQUAL(geom.ShapeFunctionsValues().size1(), 1);
            KRATOS_CHECK_EQUAL(geom.ShapeFunctionsValues().size2(), 3);

            KRATOS_CHECK_MATRIX_NEAR(
                geom.ShapeFunctionsValues(),
                p_this_quadrature_point->ShapeFunctionsValues(),
                1e-6);
        }
    } // namespace Testing
}  // namespace Kratos.
