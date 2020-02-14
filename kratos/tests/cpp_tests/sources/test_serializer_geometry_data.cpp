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
#include "geometries/point.h"

#include "includes/node.h"

#include "includes/model_part.h"
#include "containers/model.h"

#include "geometries/line_2d_3.h"
#include "geometries/quadrature_point_geometry.h"
#include "geometries/triangle_2d_3.h"

namespace Kratos {
    namespace Testing {
        typedef Geometry<Node<3>> GeometryType;
        typedef Node<3> NodeType;

        Geometry<Node<3>>::Pointer GenerateLine2D3Test(ModelPart& rModelPart) {
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.5, 0.0, 0.0);

            return Kratos::make_shared<Line2D3<Node<3>>>(
                rModelPart.pGetNode(1),
                rModelPart.pGetNode(2),
                rModelPart.pGetNode(3)
                );
        }

        Geometry<Node<3>>::Pointer GenerateLine2D3TestLoad(ModelPart& rModelPart) {
            rModelPart.CreateNewNode(4, 20, 0.0, 0.0);
            rModelPart.CreateNewNode(5, 20, 0.0, 0.0);
            rModelPart.CreateNewNode(6, 20, 0.0, 0.0);

            return Kratos::make_shared<Line2D3<Node<3>>>(
                rModelPart.pGetNode(4),
                rModelPart.pGetNode(5),
                rModelPart.pGetNode(6)
                );
        }

        /** Generates a point type sample triangle2D3.
         * Generates a point type right triangle with origin in the origin and leg size 1.
         * @return  Pointer to a triangle2D3
         */
        Geometry<Node<3>>::Pointer GeneratePointsTriangle2D3TestQP(ModelPart& rModelPart) {
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

            return Kratos::make_shared<Triangle2D3<NodeType>>(
                rModelPart.pGetNode(1),
                rModelPart.pGetNode(2),
                rModelPart.pGetNode(3)
                );
        }

        Geometry<Node<3>>::Pointer GenerateQuadraturePoint2(ModelPart& rModelPart) {
            auto triangle = GeneratePointsTriangle2D3TestQP(rModelPart);

            auto integration_points = triangle->IntegrationPoints();

            auto r_N = triangle->ShapeFunctionsValues();

            Matrix N_i = ZeroMatrix(1, triangle->size());
            for (IndexType j = 0; j < triangle->size(); ++j)
            {
                N_i(0, j) = r_N(0, j);
            }
            Matrix DN_De = triangle->ShapeFunctionLocalGradient(0);

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                GeometryData::GI_GAUSS_1,
                integration_points[0],
                N_i,
                DN_De);

            Geometry<Node<3>>::Pointer p_this_quadrature_point(
                Kratos::make_shared<QuadraturePointGeometry<Node<3>, 2, 2>>(
                    triangle->Points(),
                    data_container,
                    triangle.get()));

            return p_this_quadrature_point;
        }

        KRATOS_TEST_CASE_IN_SUITE(SerializerLine2D3, KratosCoreFastSuite)
        {
            Model model;
            auto& mp = model.CreateModelPart("SerializerLine2D3");

            StreamSerializer serializer;

            auto line_saved = GenerateLine2D3Test(mp);

            auto line_loaded = GenerateLine2D3TestLoad(mp);

            const std::string tag_string("LoadLine");

            KRATOS_CHECK_EQUAL(line_saved->size(), 3);
            serializer.save(tag_string, line_saved);
            serializer.load(tag_string, line_loaded);

            KRATOS_CHECK_NEAR((*line_saved)[0].X(), (*line_loaded)[0].X(), 1e-6);
            KRATOS_CHECK_NEAR((*line_saved)[0].Y(), (*line_loaded)[0].Y(), 1e-6);
            KRATOS_CHECK_NEAR((*line_saved)[0].Z(), (*line_loaded)[0].Z(), 1e-6);

            KRATOS_CHECK_MATRIX_NEAR(line_saved->ShapeFunctionsValues(), line_loaded->ShapeFunctionsValues(), 1e-6);

            KRATOS_CHECK_EQUAL(line_saved->size(), line_loaded->size());
            KRATOS_CHECK_EQUAL(line_saved->GetDefaultIntegrationMethod(), line_loaded->GetDefaultIntegrationMethod());

            KRATOS_CHECK_EQUAL(&(line_saved->GetGeometryData()), &(line_loaded->GetGeometryData()));
        }
    } // namespace Testing
}  // namespace Kratos.
