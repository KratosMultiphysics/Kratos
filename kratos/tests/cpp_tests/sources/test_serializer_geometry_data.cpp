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
        typedef Geometry<Node> GeometryType;
        typedef Node NodeType;

        Geometry<Node>::Pointer GenerateLine2D3Test(ModelPart& rModelPart) {
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.5, 0.0, 0.0);

            return Kratos::make_shared<Line2D3<Node>>(
                rModelPart.pGetNode(1),
                rModelPart.pGetNode(2),
                rModelPart.pGetNode(3)
                );
        }

        Geometry<Node>::Pointer GenerateLine2D3TestLoad(ModelPart& rModelPart) {
            rModelPart.CreateNewNode(4, 20, 0.0, 0.0);
            rModelPart.CreateNewNode(5, 20, 0.0, 0.0);
            rModelPart.CreateNewNode(6, 20, 0.0, 0.0);

            return Kratos::make_shared<Line2D3<Node>>(
                rModelPart.pGetNode(4),
                rModelPart.pGetNode(5),
                rModelPart.pGetNode(6)
                );
        }

        /** Generates a point type sample triangle2D3.
         * Generates a point type right triangle with origin in the origin and leg size 1.
         * @return  Pointer to a triangle2D3
         */
        Geometry<Node>::Pointer GeneratePointsTriangle2D3TestQP(ModelPart& rModelPart) {
            rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
            rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
            rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);

            return Kratos::make_shared<Triangle2D3<NodeType>>(
                rModelPart.pGetNode(1),
                rModelPart.pGetNode(2),
                rModelPart.pGetNode(3)
                );
        }

        Geometry<Node>::Pointer GenerateQuadraturePoint2DTest(ModelPart& rModelPart) {
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
                GeometryData::IntegrationMethod::GI_GAUSS_1,
                integration_points[0],
                N_i,
                DN_De);

            Geometry<Node>::Pointer p_this_quadrature_point(
                Kratos::make_shared<QuadraturePointGeometry<Node, 2, 2>>(
                    triangle->Points(),
                    data_container,
                    triangle.get()));

            return p_this_quadrature_point;
        }

        Geometry<Node>::Pointer GenerateQuadraturePoint2DTestLoad() {

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                GeometryData::IntegrationMethod::GI_GAUSS_1,
                IntegrationPoint<3>(),
                Matrix(),
                Matrix());

            Geometry<Node>::Pointer p_this_quadrature_point(
                Kratos::make_shared<QuadraturePointGeometry<Node, 2>>(
                    GeometryType::PointsArrayType(3),
                    data_container));

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

        KRATOS_TEST_CASE_IN_SUITE(SerializerQuadraturePoint, KratosCoreFastSuite)
        {
            Model model;
            auto& mp = model.CreateModelPart("SerializerQuadraturePoint");

            StreamSerializer serializer;

            auto quadrature_saved = GenerateQuadraturePoint2DTest(mp);

            auto quadrature_loaded = GenerateQuadraturePoint2DTestLoad(); // Empty quadrature point geometry

            serializer.save("qp", quadrature_saved);
            serializer.load("qp", quadrature_loaded);

            // Check coordinates
            for(unsigned int i = 0; i < 3; ++i){
                KRATOS_CHECK_NEAR((*quadrature_saved)[i].X(), (*quadrature_loaded)[i].X(), 1e-6);
                KRATOS_CHECK_NEAR((*quadrature_saved)[i].Y(), (*quadrature_loaded)[i].Y(), 1e-6);
                KRATOS_CHECK_NEAR((*quadrature_saved)[i].Z(), (*quadrature_loaded)[i].Z(), 1e-6);
            }
            // Check size and dimension
            KRATOS_CHECK_EQUAL(quadrature_saved->size(), quadrature_loaded->size());
            KRATOS_CHECK_EQUAL(quadrature_saved->WorkingSpaceDimension(), quadrature_loaded->WorkingSpaceDimension());
            KRATOS_CHECK_EQUAL(quadrature_saved->GetDefaultIntegrationMethod(), quadrature_loaded->GetDefaultIntegrationMethod());

            // Check integration point
            IntegrationPoint<3> point_loaded = quadrature_loaded->GetGeometryData().IntegrationPoints()[0];
            IntegrationPoint<3> point_saved = quadrature_saved->GetGeometryData().IntegrationPoints()[0];
            KRATOS_CHECK_NEAR(point_saved.X(), point_loaded.X(), 1e-6);
            KRATOS_CHECK_NEAR(point_saved.Y(), point_loaded.Y(), 1e-6);
            KRATOS_CHECK_NEAR(point_saved.Z(), point_loaded.Z(), 1e-6);
            KRATOS_CHECK_NEAR(point_saved.Weight(), point_loaded.Weight(),1e-6);

            // Check shape functions values
            KRATOS_CHECK_MATRIX_NEAR(quadrature_saved->ShapeFunctionsValues(),
                                     quadrature_loaded->ShapeFunctionsValues(), 1e-6);

            // Check shape functions local gradients
            KRATOS_CHECK_EQUAL(quadrature_loaded->GetGeometryData().ShapeFunctionsLocalGradients().size(),1);
            KRATOS_CHECK_MATRIX_NEAR(quadrature_saved->GetGeometryData().ShapeFunctionsLocalGradients()[0],
                                     quadrature_loaded->GetGeometryData().ShapeFunctionsLocalGradients()[0], 1e-6);

            // Check Dimensions of geometry data
            KRATOS_CHECK_EQUAL(quadrature_saved->GetGeometryData().WorkingSpaceDimension(), quadrature_loaded->GetGeometryData().WorkingSpaceDimension());
            KRATOS_CHECK_EQUAL(quadrature_saved->GetGeometryData().LocalSpaceDimension(), quadrature_loaded->GetGeometryData().LocalSpaceDimension());
        }
    } // namespace Testing
}  // namespace Kratos.
