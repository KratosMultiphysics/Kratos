//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//

// Project includes
#include "brep_face.h"

namespace Kratos
{
    void BrepFace::GetGeometryIntegration(ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        Properties::Pointer this_property = rModelPart.pGetProperties(rPropertiesId);

        auto spans_u = m_node_surface_geometry_3d->SpansU();
        auto spans_v = m_node_surface_geometry_3d->SpansV();

        ANurbs::SurfaceShapeEvaluator<double> shape(
            m_node_surface_geometry_3d->DegreeU(),
            m_node_surface_geometry_3d->DegreeV(),
            rShapeFunctionDerivativesOrder);

        for (int i = 0; i < spans_u.size(); ++i)
        {
            for (int j = 0; j < spans_v.size(); ++j)
            {
                ANurbs::Interval<double> domain_u(spans_u[i].T0(), spans_u[i].T1());
                ANurbs::Interval<double> domain_v(spans_v[j].T0(), spans_v[j].T1());


                //std::cout << "domain_u: " << spans_u[i].T0() << ", " << spans_u[i].T1() << std::endl;
                //std::cout << "domain_v: " << spans_v[j].T0() << ", " << spans_v[j].T1() << std::endl;

                auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
                    m_node_surface_geometry_3d->DegreeU() + 1,
                    m_node_surface_geometry_3d->DegreeV() + 1,
                    domain_u,
                    domain_v);
                //std::cout << "i: " << i << ", j" << j << std::endl;
                for (int k = 0; k < integration_points.size(); ++k)
                {
                    shape.Compute(
                        m_node_surface_geometry_3d->KnotsU(),
                        m_node_surface_geometry_3d->KnotsV(),
                        m_node_surface_geometry_3d->Weights(),
                        integration_points[k].u, 
                        integration_points[k].v);

                    Element::GeometryType::PointsArrayType non_zero_control_points;

                    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
                    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
                    Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

                    for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
                    {
                        int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                        int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                        non_zero_control_points.push_back(m_node_surface_geometry_3d->Node(indexU, indexV));

                        N_0[n] = shape(0, indexU, indexV);
                        N_1(n, 0) = shape(1, indexU, indexV);
                        N_1(n, 1) = shape(2, indexU, indexV);
                        N_2(n, 0) = shape(3, indexU, indexV);
                        N_2(n, 1) = shape(5, indexU, indexV);
                        N_2(n, 2) = shape(4, indexU, indexV);
                    }

                    KRATOS_WATCH(N_0)
                    KRATOS_WATCH(N_1)
                    KRATOS_WATCH(N_2)

                    if (rType == "element")
                    {
                        int id = 0;
                        if (rModelPart.GetRootModelPart().Elements().size() > 0)
                            id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                        auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, this_property);

                        element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
                        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
                        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_2);
                        element->SetValue(INTEGRATION_WEIGHT, integration_points[k].weight);
                    }
                }
            }
        }
    }

    void BrepFace::GetGeometryIntegrationTrimmed(ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
 /*       auto clipper = ANurbs::TrimmedSurfaceClipping(0.01, 0.00001);

        clipper.clear();

        for (int l = 0; l < m_embedded_loops.size(); ++l)
        {
            clipper.BeginLoop();

            auto trimming_curves = m_embedded_loops[l].GetTrimmingCurves();
            for (int t = 0; t < trimming_curves.size(); ++t)
            {
                clipper.AddCurve(*trimming_curves[t].GetCurve2D());
            }

            clipper.EndLoop();
        }

        clipper.Compute(m_node_surface_geometry_3d);

        int degree_u = m_node_surface_geometry_3d->DegreeU();
        int degree_v = m_node_surface_geometry_3d->DegreeV();

        int degree = std::max(degree_u, degree_v) + 1;

        for (int i = 0; i < clipper.NbSpansU().size(); ++i)
        {
            for (int j = 0; j < clipper.NbSpansU().size(); ++j)
            {
                if (clipper.SpanTrimType(i, j) == ANurbs::Empty)
                    continue;
                if (clipper.SpanTrimType(i, j) == ANurbs::TrimTypes->Full)
                {
                    auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
                        degree_u + 1,
                        degree_v + 1,
                        clipper.SpanU(i),
                        clipper.SpanV(j));

                    CreateIntegrationElementsConditions(integration_points, m_model_part, rType, rName,
                        rPropertiesId, rShapeFunctionDerivativesOrder, rVariables);

                    continue;
                }
                else
                {
                    for (int p = 0; p < clipper.SpanPolygons(i, j); ++p)
                    {
                        auto integration_point_polygon = ANurbs::PolygonIntegrationPoints();

                        integration_point_polygon.Compute(degree, clipper.SpanPolygons(i, j)[p]);

                        std::vector<ANurbs::IntegrationPoint2<double>> integration_points =
                            integration_point_polygon.IntegrationPoints();

                        CreateIntegrationElementsConditions(integration_points, m_model_part, rType, rName,
                            rPropertiesId, rShapeFunctionDerivativesOrder, rVariables);
                    }
                }
            }
        }*/

    }

    void BrepFace::CreateIntegrationElementsConditions(
        std::vector<ANurbs::IntegrationPoint2<double>> rIntegrationPoints,
        ModelPart& rModelPart,
        const std::string& rType,
        const std::string& rName,
        const int& rPropertiesId,
        const int& rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        Properties::Pointer this_property = rModelPart.pGetProperties(rPropertiesId);

        ANurbs::SurfaceShapeEvaluator<double> shape(
            m_node_surface_geometry_3d->DegreeU(),
            m_node_surface_geometry_3d->DegreeV(),
            rShapeFunctionDerivativesOrder);

        for (int k = 0; k < rIntegrationPoints.size(); ++k)
        {
            shape.Compute(
                m_node_surface_geometry_3d->KnotsU(),
                m_node_surface_geometry_3d->KnotsV(),
                m_node_surface_geometry_3d->Weights(),
                rIntegrationPoints[k].u,
                rIntegrationPoints[k].v);

            Element::GeometryType::PointsArrayType non_zero_control_points;

            Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
            Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
            Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

            for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
            {
                int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                non_zero_control_points.push_back(m_node_surface_geometry_3d->Node(indexU, indexV));

                N_0[n] = shape(0, indexU, indexV);
                N_1(n, 0) = shape(1, indexU, indexV);
                N_1(n, 1) = shape(2, indexU, indexV);
                N_2(n, 0) = shape(3, indexU, indexV);
                N_2(n, 1) = shape(5, indexU, indexV);
                N_2(n, 2) = shape(4, indexU, indexV);
            }

            KRATOS_WATCH(N_0)
            KRATOS_WATCH(N_1)
            KRATOS_WATCH(N_2)

            if (rType == "element")
            {
                int id = 0;
                if (rModelPart.GetRootModelPart().Elements().size() > 0)
                    id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, this_property);

                element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
                element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
                element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_2);
                element->SetValue(INTEGRATION_WEIGHT, rIntegrationPoints[k].weight);
            }
        }
    }


    ///Constructor
    BrepFace::BrepFace(
        int rBrepId,
        bool rIsTrimmed,
        bool rIsRational,
        std::vector<BrepBoundaryLoop>& rTrimmingLoops,
        std::vector<BrepBoundaryLoop>& rEmbeddedLoops,
        std::vector<EmbeddedPoint>& rEmbeddedPoints,
        Vector& rKnotVectorU,
        Vector& rKnotVectorV,
        int& rP,
        int& rQ,
        IntVector& rControlPointIds,
        ModelPart& rModelPart)
        : m_trimming_loops(rTrimmingLoops),
          m_is_trimmed(rIsTrimmed),
          m_is_rational(rIsRational),
          m_embedded_loops(rEmbeddedLoops),
          m_embedded_points(rEmbeddedPoints),
          m_control_points_ids(rControlPointIds),
          m_model_part(rModelPart),
          IndexedObject(rBrepId),
          Flags()
    {
        int number_of_nodes_u = rKnotVectorU.size() - rP - 1;
        int number_of_nodes_v = rKnotVectorV.size() - rQ - 1;

        m_node_surface_geometry_3d = Kratos::make_unique<NodeSurfaceGeometry3D>(
            rP, rQ, number_of_nodes_u, number_of_nodes_v);

        for (int i = 0; i < rKnotVectorU.size()-2; ++i)
        {
            m_node_surface_geometry_3d->SetKnotU(i, rKnotVectorU(i+1));
        }

        for (int i = 0; i < rKnotVectorV.size()-2; ++i)
        {
            m_node_surface_geometry_3d->SetKnotV(i, rKnotVectorV(i+1));
        }

        for (int i = 0; i < number_of_nodes_u; ++i)
        {
            for (int j = 0; j < number_of_nodes_v; ++j)
            {
                Node<3>::Pointer node = rModelPart.pGetNode(rControlPointIds[i*(number_of_nodes_u)+j]);
                m_node_surface_geometry_3d->SetNode(i, j, node);
                if (rIsRational)
                {
                    m_node_surface_geometry_3d->SetWeight(i, j, node->GetValue(NURBS_CONTROL_POINT_WEIGHT));
                }
            }
        }

        //KRATOS_WATCH(m_node_surface_geometry_3d->KnotsU())
        //KRATOS_WATCH(m_node_surface_geometry_3d->KnotsV())
        //KRATOS_WATCH(m_node_surface_geometry_3d->NbPoles())
        //KRATOS_WATCH(m_node_surface_geometry_3d->Weights())
    }
} // namespace Kratos.