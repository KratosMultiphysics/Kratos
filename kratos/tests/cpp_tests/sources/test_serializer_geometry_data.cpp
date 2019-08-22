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

#include "includes/node.h"

#include "geometries/line_2d_3.h"
#include "geometries/quadrature_point.h"
#include "geometries/triangle_2d_3.h"

namespace Kratos {
    namespace Testing {
        typedef Geometry<Node<3>> GeometryType;
        typedef Node<3> NodeType;

        Geometry<Node<3>>::Pointer GenerateLine2D3() {
            return Kratos::make_shared<Line2D3<Node<3>>>(
                new NodeType(1, 0.0, 0.0, 0.0),
                new NodeType(2, 1.0, 0.0, 0.0),
                new NodeType(3, 0.5, 0.0, 0.0)
                );
        }

        Geometry<Node<3>>::Pointer GenerateLine2D32() {
            return Kratos::make_shared<Line2D3<Node<3>>>(
                new NodeType(4, 20, 0.0, 0.0),
                new NodeType(5, 20, 0.0, 0.0),
                new NodeType(6, 20, 0.0, 0.0)
                );
        }

        /** Generates a point type sample triangle2D3.
 * Generates a point type right triangle with origin in the origin and leg size 1.
 * @return  Pointer to a triangle2D3
 */
        Geometry<Node<3>>::Pointer GeneratePointsTriangle2D32() {
            return Kratos::make_shared<Triangle2D3<NodeType>>(
                new NodeType(1, 0.0, 0.0, 0.0),
                new NodeType(2, 1.0, 0.0, 0.0),
                new NodeType(3, 0.0, 1.0, 0.0)
                );
        }

        Geometry<Node<3>>::Pointer GenerateQuadraturePoint2() {
            auto triangle = GeneratePointsTriangle2D32();

            auto integration_points = triangle->IntegrationPoints();

            auto r_N = triangle->ShapeFunctionsValues();

                Matrix N_i = ZeroMatrix(1, triangle->size());
                for (std::size_t j = 0; j < triangle->size(); ++j)
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

                Geometry<Node<3>>::Pointer p_this_quadrature_point(
                    Kratos::make_shared<QuadraturePoint<2, 2, Node<3>>>(
                        triangle->Points(),
                        ips_container,
                        N_container,
                        DN_De_container));

                return p_this_quadrature_point;
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerLine2D3, KratosCoreFastSuite)
        {
            StreamSerializer serializer;

            auto line_saved = GenerateLine2D3();

            auto line_loaded = GenerateLine2D32();

            const std::string tag_string("LoadLine");

            serializer.save(tag_string, line_saved);
            serializer.load(tag_string, line_loaded);

            KRATOS_CHECK_NEAR((*line_saved)[0].X(), (*line_loaded)[0].X(), 1e-6);
            KRATOS_CHECK_NEAR((*line_saved)[0].Y(), (*line_loaded)[0].Y(), 1e-6);
            KRATOS_CHECK_NEAR((*line_saved)[0].Z(), (*line_loaded)[0].Z(), 1e-6);

            for (std::size_t i = 0; i < line_saved->size(); ++i)
            {
                for (std::size_t j = 0; j < line_saved->IntegrationPointsNumber(); ++j)
                {
                    KRATOS_CHECK_NEAR(line_saved->ShapeFunctionValue(j,i), line_loaded->ShapeFunctionValue(j, i), 1e-6);
                }
            }

            KRATOS_CHECK_EQUAL(line_saved->size(), line_loaded->size());
            KRATOS_CHECK_EQUAL(line_saved->GetDefaultIntegrationMethod(), line_loaded->GetDefaultIntegrationMethod());
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerQuadraturePoint, KratosCoreFastSuite)
        {
            StreamSerializer serializer;

            auto quadrature_saved = GenerateQuadraturePoint2();

            auto quadrature_loaded = GenerateQuadraturePoint2();

            const std::string tag_string("LoadQuadraturePoint");

            serializer.save(tag_string, quadrature_saved);
            serializer.load(tag_string, quadrature_loaded);

            KRATOS_CHECK_NEAR((*quadrature_saved)[0].X(), (*quadrature_loaded)[0].X(), 1e-6);
            KRATOS_CHECK_NEAR((*quadrature_saved)[0].Y(), (*quadrature_loaded)[0].Y(), 1e-6);
            KRATOS_CHECK_NEAR((*quadrature_saved)[0].Z(), (*quadrature_loaded)[0].Z(), 1e-6);

            for (std::size_t i = 0; i < quadrature_saved->size(); ++i)
            {
                for (std::size_t j = 0; j < quadrature_saved->IntegrationPointsNumber(); ++j)
                {
                    KRATOS_CHECK_NEAR(quadrature_saved->ShapeFunctionValue(j, i), quadrature_loaded->ShapeFunctionValue(j, i), 1e-6);
                }
            }

            KRATOS_CHECK_EQUAL(quadrature_saved->size(), quadrature_loaded->size());
            KRATOS_CHECK_EQUAL(quadrature_saved->WorkingSpaceDimension(), quadrature_loaded->WorkingSpaceDimension());
            KRATOS_CHECK_EQUAL(quadrature_saved->GetDefaultIntegrationMethod(), quadrature_loaded->GetDefaultIntegrationMethod());
        }
    } // namespace Testing
}  // namespace Kratos.
