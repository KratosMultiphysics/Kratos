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

                auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
                    m_node_surface_geometry_3d->DegreeU() + 1,
                    m_node_surface_geometry_3d->DegreeV() + 1,
                    domain_u,
                    domain_v);

                for (int k = 0; k < integration_points.size(); ++k)
                {
                    shape.Compute(
                        m_node_surface_geometry_3d->KnotsU(),
                        m_node_surface_geometry_3d->KnotsV(),
                        //m_node_surface_geometry_3d->Weights(),
                        integration_points[k].u, 
                        integration_points[k].v);

                    Element::GeometryType::PointsArrayType non_zero_control_points;

                    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
                    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
                    Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

                    Vector coords = ZeroVector(3);
                    KRATOS_WATCH(integration_points[k].u)
                    KRATOS_WATCH(integration_points[k].v)
                    KRATOS_WATCH(integration_points[k].weight)
                    //KRATOS_WATCH(shape.NonzeroPoleIndices().size())
                    for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
                    {
                        int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                        int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                        //std::cout << "indexU: " << indexU << ", indexV: " << indexV << std::endl;
                        //std::cout << "shape.NonzeroPoleIndices()[n].first: " << shape.NonzeroPoleIndices()[n].first 
                        //    << ", shape.NonzeroPoleIndices()[n].second: " << shape.NonzeroPoleIndices()[n].second << std::endl;

                        //std::cout << m_node_surface_geometry_3d->Node(shape.NonzeroPoleIndices()[n].first,
                        //    shape.NonzeroPoleIndices()[n].second)->Id() << std::endl;

                        non_zero_control_points.push_back(m_node_surface_geometry_3d->Node(
                            shape.NonzeroPoleIndices()[n].first, shape.NonzeroPoleIndices()[n].second));

                        N_0[n] = shape(0, indexU, indexV);
                        N_1(n, 0) = shape(1, indexU, indexV);
                        N_1(n, 1) = shape(2, indexU, indexV);
                        N_2(n, 0) = shape(3, indexU, indexV);
                        N_2(n, 1) = shape(5, indexU, indexV);
                        N_2(n, 2) = shape(4, indexU, indexV);

                        coords[0] += N_0[n] * m_node_surface_geometry_3d->Node(
                            shape.NonzeroPoleIndices()[n].first, shape.NonzeroPoleIndices()[n].second)->X();
                        coords[1] += N_0[n] * m_node_surface_geometry_3d->Node(
                            shape.NonzeroPoleIndices()[n].first, shape.NonzeroPoleIndices()[n].second)->Y();
                        coords[2] += N_0[n] * m_node_surface_geometry_3d->Node(
                            shape.NonzeroPoleIndices()[n].first, shape.NonzeroPoleIndices()[n].second)->Z();
                    }
                    //KRATOS_WATCH(coords)
                    KRATOS_WATCH(N_0)
                    //KRATOS_WATCH(N_1)
                    //KRATOS_WATCH(N_2)

                    if (rType == "element")
                    {
                        int id = 0;
                        if (rModelPart.GetRootModelPart().Elements().size() > 0)
                            id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                        rModelPart.AddNodes(non_zero_control_points.begin(), non_zero_control_points.end());

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
        //auto clipper = ANurbs::TrimmedSurfaceClipping<Kratos::array_1d<double, 2>>(0.01, 0.00001);

        //clipper.Clear();

        //for (int l = 0; l < m_embedded_loops.size(); ++l)
        //{
        //    clipper.BeginLoop();

        //    auto trimming_curves = m_embedded_loops[l].GetTrimmingCurves();
        //    for (int t = 0; t < trimming_curves.size(); ++t)
        //    {
        //        auto crv = (trimming_curves[t].GetCurve2DNonconst());
        //        //clipper.AddCurve(*crv);
        //    }

        //    clipper.EndLoop();
        //}

        //clipper.Compute(m_node_surface_geometry_3d->SpansU(), m_node_surface_geometry_3d->SpansV());

        //int degree_u = m_node_surface_geometry_3d->DegreeU();
        //int degree_v = m_node_surface_geometry_3d->DegreeV();

        //int degree = std::max(degree_u, degree_v) + 1;

        //for (int i = 0; i < clipper.NbSpansU(); ++i)
        //{
        //    for (int j = 0; j < clipper.NbSpansU(); ++j)
        //    {
        //        if (clipper.SpanTrimType(i, j) == ANurbs::Empty)
        //            continue;
        //        if (clipper.SpanTrimType(i, j) == ANurbs::Full)
        //        {
        //            auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
        //                degree_u + 1,
        //                degree_v + 1,
        //                clipper.SpanU(i),
        //                clipper.SpanV(j));

        //            CreateIntegrationElementsConditions(integration_points, m_model_part, rType, rName,
        //                rPropertiesId, rShapeFunctionDerivativesOrder, rVariables);

        //            continue;
        //        }
        //        else
        //        {
        //            auto polygons = clipper.SpanPolygons(i, j);
        //            for (int p = 0; p < polygons.size(); ++p)
        //            {
        //                auto integration_point_polygon = ANurbs::PolygonIntegrationPoints<Kratos::array_1d<double, 2>>();

        //                integration_point_polygon.Compute(degree, polygons[p]);

        //                std::vector<ANurbs::IntegrationPoint2<double>> integration_points;
        //                for (int integration_point = 0; integration_point < integration_point_polygon.NbIntegrationPoints(); ++integration_point)
        //                {
        //                    integration_points[integration_point] = integration_point_polygon.IntegrationPoint(integration_point);
        //                }
        //                CreateIntegrationElementsConditions(integration_points, m_model_part, rType, rName,
        //                    rPropertiesId, rShapeFunctionDerivativesOrder, rVariables);
        //            }
        //        }
        //    }
        //}
    }

    void BrepFace::GetIntegrationBrepEdge(
        ModelPart& rModelPart,
        const int& trim_index,
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

        auto curve_2d = GetTrimCurve(trim_index);

        //auto curveGeometry2d = Kratos::make_shared<CurveGeometry<2>>(
        //    2,      // DegreeU
        //    4,      // NumberOfPoles
        //    false   // IsRational
        //    );



        //// NURBS-Surface-Geometry based on FE-Nodes

        //auto nodeSurfaceGeometry3d = Kratos::make_shared<NodeSurfaceGeometry3D>(
        //    2,      // DegreeU
        //    2,      // DegreeV
        //    4,      // NumberOfPolesU
        //    4       // NumberOfPolesV
        //    );

        //auto curveOnSurface3d = Kratos::make_shared<CurveOnSurface<3>>(
        //    curveGeometry2d,            // CurveGeometry2D
        //    nodeSurfaceGeometry3d,      // SurfaceGeometry3D
        //    curveGeometry2d->Domain()   // Domain
        //    );

        //auto curve_on_surface_3d = Kratos::make_shared<CurveOnSurface<3>>(
        //        curve_2d->CurveGeometry(), m_node_surface_geometry_3d, curve_2d->Domain());
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

            //KRATOS_WATCH(N_0)
            //KRATOS_WATCH(N_1)
            //KRATOS_WATCH(N_2)

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

    //Kratos::Element::Pointer BrepFace::CreateIntegrationElement(
    //    const int& rU,
    //    const int& rV,
    //    ModelPart& rModelPart,
    //    const std::string& rType,
    //    const std::string& rName,
    //    Properties::Pointer pProperty,
    //    const int& rShapeFunctionDerivativesOrder,
    //    std::vector<std::string> rVariables) const
    //{
    //    ANurbs::SurfaceShapeEvaluator<double> shape(
    //        m_node_surface_geometry_3d->DegreeU(),
    //        m_node_surface_geometry_3d->DegreeV(),
    //        rShapeFunctionDerivativesOrder);

    //    shape.Compute(
    //        m_node_surface_geometry_3d->KnotsU(),
    //        m_node_surface_geometry_3d->KnotsV(),
    //        m_node_surface_geometry_3d->Weights(),
    //        rU,
    //        rV);

    //    Element::GeometryType::PointsArrayType non_zero_control_points;

    //    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
    //    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
    //    Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

    //    for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
    //    {
    //        int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
    //        int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

    //        non_zero_control_points.push_back(m_node_surface_geometry_3d->Node(indexU, indexV));

    //        N_0[n] = shape(0, indexU, indexV);
    //        N_1(n, 0) = shape(1, indexU, indexV);
    //        N_1(n, 1) = shape(2, indexU, indexV);
    //        N_2(n, 0) = shape(3, indexU, indexV);
    //        N_2(n, 1) = shape(5, indexU, indexV);
    //        N_2(n, 2) = shape(4, indexU, indexV);
    //    }

    //    if (rType == "element")
    //    {
    //        int id = 0;
    //        if (rModelPart.GetRootModelPart().Elements().size() > 0)
    //            id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

    //        auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, pProperty);

    //        element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
    //        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
    //        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_2);
    //        return element;
    //    }
    //}

    //void BrepFace::CreateIntegrationElementSlave(
    //    const int& rU,
    //    const int& rV,
    //    Kratos::Element::Pointer pElementPoint,
    //    const int& rShapeFunctionDerivativesOrder,
    //    std::vector<std::string> rVariables) const
    //{
    //    ANurbs::SurfaceShapeEvaluator<double> shape(
    //        m_node_surface_geometry_3d->DegreeU(),
    //        m_node_surface_geometry_3d->DegreeV(),
    //        rShapeFunctionDerivativesOrder);

    //    shape.Compute(
    //        m_node_surface_geometry_3d->KnotsU(),
    //        m_node_surface_geometry_3d->KnotsV(),
    //        m_node_surface_geometry_3d->Weights(),
    //        rU,
    //        rV);

    //    Element::GeometryType::PointsArrayType non_zero_control_points;

    //    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
    //    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
    //    Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

    //    for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
    //    {
    //        int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
    //        int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

    //        non_zero_control_points.push_back(m_node_surface_geometry_3d->Node(indexU, indexV));

    //        N_0[n]    = shape(0, indexU, indexV);
    //        N_1(n, 0) = shape(1, indexU, indexV);
    //        N_1(n, 1) = shape(2, indexU, indexV);
    //        N_2(n, 0) = shape(3, indexU, indexV);
    //        N_2(n, 1) = shape(5, indexU, indexV);
    //        N_2(n, 2) = shape(4, indexU, indexV);
    //    }

    //    pElementPoint->SetValue(SHAPE_FUNCTION_VALUES_SLAVE, N_0);
    //    pElementPoint->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, N_1);
    //    pElementPoint->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, N_2);
    //}

    //void BrepFace::CreateIntegrationElementSlave(
    //    const int& rU,
    //    const int& rV,
    //    Kratos::Element::Pointer pElementPoint,
    //    const int& rShapeFunctionDerivativesOrder,
    //    std::vector<std::string> rVariables) const
    //{
    //    ANurbs::SurfaceShapeEvaluator<double> shape(
    //        m_node_surface_geometry_3d->DegreeU(),
    //        m_node_surface_geometry_3d->DegreeV(),
    //        rShapeFunctionDerivativesOrder);

    //    shape.Compute(
    //        m_node_surface_geometry_3d->KnotsU(),
    //        m_node_surface_geometry_3d->KnotsV(),
    //        m_node_surface_geometry_3d->Weights(),
    //        rU,
    //        rV);

    //    Element::GeometryType::PointsArrayType non_zero_control_points;

    //    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
    //    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
    //    Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

    //    for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
    //    {
    //        int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
    //        int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

    //        non_zero_control_points.push_back(m_node_surface_geometry_3d->Node(indexU, indexV));

    //        N_0[n] = shape(0, indexU, indexV);
    //        N_1(n, 0) = shape(1, indexU, indexV);
    //        N_1(n, 1) = shape(2, indexU, indexV);
    //        N_2(n, 0) = shape(3, indexU, indexV);
    //        N_2(n, 1) = shape(5, indexU, indexV);
    //        N_2(n, 2) = shape(4, indexU, indexV);
    //    }

    //    pElementPoint->SetValue(SHAPE_FUNCTION_VALUES_SLAVE, N_0);
    //    pElementPoint->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, N_1);
    //    pElementPoint->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, N_2);
    //}

    void BrepFace::EvaluatePoint(
        const int& rU,
        const int& rV,
        Element::GeometryType::PointsArrayType& rNonZeroControlPoints,
        Vector& rShapeFunction,
        Matrix& rShapeFunctionDerivative,
        Matrix& rShapeFunctionSecondDerivative
    ) const
    {
        ANurbs::SurfaceShapeEvaluator<double> shape(
            m_node_surface_geometry_3d->DegreeU(),
            m_node_surface_geometry_3d->DegreeV(),
            3);

        shape.Compute(
            m_node_surface_geometry_3d->KnotsU(),
            m_node_surface_geometry_3d->KnotsV(),
            m_node_surface_geometry_3d->Weights(),
            rU,
            rV);

        if (rShapeFunction.size() != shape.NbNonzeroPoles())
            rShapeFunction.resize(shape.NbNonzeroPoles());
        rShapeFunction = ZeroVector(shape.NbNonzeroPoles());

        if (rShapeFunctionDerivative.size1() != shape.NbNonzeroPoles() && rShapeFunctionDerivative.size2() != 2)
            rShapeFunctionDerivative.resize(rShapeFunctionDerivative.size1(), 2);
        rShapeFunctionDerivative = ZeroMatrix(shape.NbNonzeroPoles(), 2);

        if (rShapeFunctionSecondDerivative.size1() != shape.NbNonzeroPoles() && rShapeFunctionSecondDerivative.size2() != 2)
            rShapeFunctionSecondDerivative.resize(rShapeFunctionSecondDerivative.size1(), 3);
        rShapeFunctionSecondDerivative = ZeroMatrix(shape.NbNonzeroPoles(), 3);

        for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
        {
            int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
            int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

            rNonZeroControlPoints.push_back(m_node_surface_geometry_3d->Node(indexU, indexV));

            rShapeFunction[n] = shape(0, indexU, indexV);
            rShapeFunctionDerivative(n, 0) = shape(1, indexU, indexV);
            rShapeFunctionDerivative(n, 1) = shape(2, indexU, indexV);
            rShapeFunctionSecondDerivative(n, 0) = shape(3, indexU, indexV);
            rShapeFunctionSecondDerivative(n, 1) = shape(5, indexU, indexV);
            rShapeFunctionSecondDerivative(n, 2) = shape(4, indexU, indexV);
        }
    }

    const Kratos::shared_ptr<NodeSurfaceGeometry3D> BrepFace::GetSurface() const
    {
        return m_node_surface_geometry_3d;
    }

    const Kratos::shared_ptr<ANurbs::Curve<Kratos::array_1d<double, 2>>> BrepFace::GetTrimCurve(const int& trim_index) const
    {
        for (int i = 0; i < m_trimming_loops.size(); ++i)
        {
            std::vector<BrepTrimmingCurve> trimming_loops = m_trimming_loops[i].GetTrimmingCurves();
            for (int j = 0; j < trimming_loops.size(); ++j)
            {
                if (trimming_loops[j].GetTrimIndex() == trim_index)
                    return trimming_loops[j].GetCurve2D();
            }
        }
        for (int i = 0; i < m_embedded_loops.size(); ++i)
        {
            std::vector<BrepTrimmingCurve> trimming_loops = m_trimming_loops[i].GetTrimmingCurves();
            for (int j = 0; j < trimming_loops.size(); ++j)
            {
                if (trimming_loops[j].GetTrimIndex() == trim_index)
                    return trimming_loops[j].GetCurve2D();
            }
        }
        KRATOS_ERROR << "Trimming curve of Id: " << trim_index << " does not exist." << std::endl;
    }

    //Kratos::shared_ptr<ANurbs::Curve<Kratos::array_1d<double, 2>>> BrepFace::GetTrimCurveNonConst(const int& trim_index)
    //{
    //    for (int i = 0; i < m_trimming_loops.size(); ++i)
    //    {
    //        std::vector<BrepTrimmingCurve> trimming_loops = m_trimming_loops[i].GetTrimmingCurves();
    //        for (int j = 0; j < trimming_loops.size(); ++j)
    //        {
    //            if (trimming_loops[j].GetTrimIndex() == trim_index)
    //                return trimming_loops[j].GetCurve2D();
    //        }
    //    }
    //    for (int i = 0; i < m_embedded_loops.size(); ++i)
    //    {
    //        std::vector<BrepTrimmingCurve> trimming_loops = m_trimming_loops[i].GetTrimmingCurves();
    //        for (int j = 0; j < trimming_loops.size(); ++j)
    //        {
    //            if (trimming_loops[j].GetTrimIndex() == trim_index)
    //                return trimming_loops[j].GetCurve2D();
    //        }
    //    }
    //    KRATOS_ERROR << "Trimming curve of Id: " << trim_index << " does not exist." << std::endl;
    //}

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

        //KRATOS_WATCH(number_of_nodes_u)
        //KRATOS_WATCH(number_of_nodes_v)

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
                Node<3>::Pointer node = rModelPart.pGetNode(rControlPointIds[i*(number_of_nodes_v)+j]);
                //std::cout << "index: " << i * (number_of_nodes_u)+j << std::endl;
                m_node_surface_geometry_3d->SetNode(i, j, node);
                if (rIsRational)
                {
                    //KRATOS_WATCH(node->GetValue(NURBS_CONTROL_POINT_WEIGHT))
                    m_node_surface_geometry_3d->SetWeight(i, j, node->GetValue(NURBS_CONTROL_POINT_WEIGHT));
                }
            }
        }

        KRATOS_WATCH(m_node_surface_geometry_3d->KnotsU())
        KRATOS_WATCH(m_node_surface_geometry_3d->KnotsV())
        KRATOS_WATCH(m_node_surface_geometry_3d->NbPoles())
        //KRATOS_WATCH(m_node_surface_geometry_3d->Weights())
    }
} // namespace Kratos.