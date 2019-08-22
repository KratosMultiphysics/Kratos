//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// Project includes
#include "testing/testing.h"
#include "includes/serializer.h"
#include "includes/file_serializer.h"
#include "includes/stream_serializer.h"
#include "geometries/point.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/quadrature_point.h"
#include "geometries/geometry_shape_function_container.h"

namespace Kratos {
    namespace Testing {
        typedef Node<3> NodeType;

        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        /** Generates a point type sample triangle2D3.
         * Generates a point type right triangle with origin in the origin and leg size 1.
         * @return  Pointer to a triangle2D3
         */
        Geometry<Node<3>>::Pointer GeneratePointsTriangle2D3() {
            return Kratos::make_shared<Triangle2D3<NodeType>>(
                new NodeType(1, 0.0, 0.0, 0.0),
                new NodeType(2, 1.0, 0.0, 0.0),
                new NodeType(3, 0.0, 1.0, 0.0)
                );
        }

        Geometry<Node<3>>::Pointer GenerateQuadraturePoint() {
            auto triangle = GeneratePointsTriangle2D3();

            //auto integration_method = triangle->GetDefaultIntegrationMethod();
            auto integration_points = triangle->IntegrationPoints();

            auto r_N = triangle->ShapeFunctionsValues();

            Matrix N_i = ZeroMatrix(1, triangle->size());
            for (IndexType j = 0; j < triangle->size(); ++j)
            {
                N_i(0, j) = r_N(0, j);
            }
            Matrix DN_De = triangle->ShapeFunctionLocalGradient(0);

            QuadraturePoint<2, 2, Node<3>>::IntegrationPointsArrayType ips(1);
            ips[0] = integration_points[0];

            QuadraturePoint<2, 2, Node<3>>::IntegrationPointsContainerType ips_container =
            { { ips } };

            QuadraturePoint<2, 2, Node<3>>::ShapeFunctionsValuesContainerType N_container =
            { { N_i } };

            auto DN_De_array = DenseVector<Matrix>(1);
            DN_De_array[0] = DN_De;

            QuadraturePoint<2, 2, Node<3>>::ShapeFunctionsLocalGradientsContainerType DN_De_container =
            { { DN_De_array } };

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                GeometryData::GI_GAUSS_1,
                integration_points[0],
                N_i,
                DN_De);

            Geometry<Node<3>>::Pointer p_this_quadrature_point(
                Kratos::make_shared<QuadraturePoint<2, 2, Node<3>>>(
                    triangle->Points(),
                    data_container,
                    triangle.get()));

            return p_this_quadrature_point;
        }

        KRATOS_TEST_CASE_IN_SUITE(QuadraturePoint2d, KratosCoreFastSuite)
        {
            auto p_this_quadrature_point = GenerateQuadraturePoint();

            KRATOS_CHECK_EQUAL(p_this_quadrature_point->size(), 3);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->WorkingSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->LocalSpaceDimension(), 2);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->Dimension(), 2);

            KRATOS_CHECK_EQUAL(p_this_quadrature_point->ShapeFunctionsValues().size1(), 1);
            KRATOS_CHECK_EQUAL(p_this_quadrature_point->ShapeFunctionsValues().size2(), 3);

            KRATOS_WATCH(p_this_quadrature_point->ShapeFunctionsValues())
            KRATOS_WATCH(p_this_quadrature_point->ShapeFunctionsLocalGradients()[0])
            KRATOS_WATCH(p_this_quadrature_point->Center())
        }

    } // namespace Testing
}  // namespace Kratos.
