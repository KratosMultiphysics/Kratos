//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//               Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//

// Project includes
#include "brep_face.h"

namespace Kratos
{
    void BrepFace::GetGeometryNodes(
        ModelPart& rModelPart,
        const int& rU,
        const int& rV) const
    {
        int number_of_cps_u = mNodeSurfaceGeometry3D->NbPolesU();
        int number_of_cps_v = mNodeSurfaceGeometry3D->NbPolesV();

        int u_start = 0;
        int u_end = number_of_cps_u;
        int v_start = 0;
        int v_end = number_of_cps_v;

        if (rU >= 0)
        {
            u_start = rU * (number_of_cps_u - 1);
            u_end = rU * (number_of_cps_u - 1) + 1;
        }
        if (rV >= 0)
        {
            v_start = rV * (number_of_cps_v - 1);
            v_end = rV * (number_of_cps_v - 1) + 1;
        }

        for (int i = u_start; i < u_end; ++i)
        {
            for (int j = v_start; j < v_end; ++j)
            {
                rModelPart.AddNode(mNodeSurfaceGeometry3D->GetNode(i, j));
            }
        }
    }

    void BrepFace::GetGeometryVariationNodes(
        ModelPart& rModelPart,
        const int& rU,
        const int& rV) const
    {
        int number_of_cps_u = mNodeSurfaceGeometry3D->NbPolesU();
        int number_of_cps_v = mNodeSurfaceGeometry3D->NbPolesV();

        int u_start = 0;
        int u_end = number_of_cps_u;
        int v_start = 0;
        int v_end = number_of_cps_v;

        if (rU == 0)
        {
            u_start = 1;
            u_end = 2;
        }
        if (rU == 1)
        {
            u_start = number_of_cps_u - 2;
            u_end = number_of_cps_u - 1;
        }
        if (rV == 0)
        {
            v_start = 1;
            v_end = 2;
        }
        if (rV == 1)
        {
            v_start = number_of_cps_v - 2;
            v_end = number_of_cps_v - 1;
        }

        for (int i = u_start; i < u_end; ++i)
        {
            for (int j = v_start; j < v_end; ++j)
            {
                rModelPart.AddNode(mNodeSurfaceGeometry3D->GetNode(i, j));
            }
        }
    }

    //void BrepFace::GetGeometryIntegration(ModelPart& rModelPart,
    //    const std::string& rType,
    //    const std::string& rName,
    //    const int rShapeFunctionDerivativesOrder,
    //    std::vector<std::string> rVariables)
    //{
    //    auto spans_u = mNodeSurfaceGeometry3D->SpansU();
    //    auto spans_v = mNodeSurfaceGeometry3D->SpansV();

    //    ANurbs::SurfaceShapeEvaluator<double> shape(
    //        mNodeSurfaceGeometry3D->DegreeU(),
    //        mNodeSurfaceGeometry3D->DegreeV(),
    //        rShapeFunctionDerivativesOrder);


    //    for (int i = 0; i < spans_u.size(); ++i)
    //    {
    //        for (int j = 0; j < spans_v.size(); ++j)
    //        {
    //            ANurbs::Interval<double> domain_u(spans_u[i].T0(), spans_u[i].T1());
    //            ANurbs::Interval<double> domain_v(spans_v[j].T0(), spans_v[j].T1());

    //            auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
    //                mNodeSurfaceGeometry3D->DegreeU() + 1,
    //                mNodeSurfaceGeometry3D->DegreeV() + 1,
    //                domain_u,
    //                domain_v);

    //            for (int k = 0; k < integration_points.size(); ++k)
    //            {
    //                array_1d<double, 2> local_coordinates;
    //                local_coordinates[0] = integration_points[k].u;
    //                local_coordinates[1] = integration_points[k].v;

    //                shape.Compute(
    //                    mNodeSurfaceGeometry3D->KnotsU(),
    //                    mNodeSurfaceGeometry3D->KnotsV(),
    //                    //mNodeSurfaceGeometry3D->Weights(),
    //                    integration_points[k].u,
    //                    integration_points[k].v);

    //                Element::GeometryType::PointsArrayType non_zero_control_points;

    //                Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
    //                Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
    //                Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

    //                Vector coords = ZeroVector(3);
    //                for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
    //                {
    //                    int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
    //                    int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

    //                    non_zero_control_points.push_back(mNodeSurfaceGeometry3D->GetNode(
    //                        shape.NonzeroPoleIndices()[n].first, shape.NonzeroPoleIndices()[n].second));

    //                    N_0[n] = shape(0, indexU, indexV);
    //                    N_1(n, 0) = shape(1, indexU, indexV);
    //                    N_1(n, 1) = shape(2, indexU, indexV);
    //                    N_2(n, 0) = shape(3, indexU, indexV);
    //                    N_2(n, 1) = shape(5, indexU, indexV);
    //                    N_2(n, 2) = shape(4, indexU, indexV);
    //                }

    //                if (rType == "element")
    //                {
    //                    int id = 1;
    //                    if (rModelPart.GetRootModelPart().Elements().size() > 0)
    //                        id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

    //                    rModelPart.AddNodes(non_zero_control_points.begin(), non_zero_control_points.end());

    //                    auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, 0);

    //                    element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
    //                    element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
    //                    element->SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, N_2);
    //                    element->SetValue(INTEGRATION_WEIGHT, integration_points[k].weight);

    //                    element->SetValue(LOCAL_COORDINATES, local_coordinates);
    //                    element->SetValue(BREP_ID, this->Id());
    //                }

    //                if (rType == "condition")
    //                {
    //                    int id = 1;
    //                    if (rModelPart.GetRootModelPart().Conditions().size() > 0)
    //                        id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;

    //                    rModelPart.AddNodes(non_zero_control_points.begin(), non_zero_control_points.end());

    //                    auto condition = rModelPart.CreateNewCondition(rName, id, non_zero_control_points, 0);

    //                    condition->SetValue(SHAPE_FUNCTION_VALUES, N_0);
    //                    condition->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
    //                    condition->SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, N_2);
    //                    condition->SetValue(INTEGRATION_WEIGHT, integration_points[k].weight);

    //                    condition->SetValue(LOCAL_COORDINATES, local_coordinates);
    //                    condition->SetValue(BREP_ID, this->Id());
    //                }
    //            }
    //        }
    //    }
    //}

    TrimmedSurfaceClipping BrepFace::GetSurfaceClipper(
        const double& rAccuracy,
        const double& rUnit)
    {
        KRATOS_WATCH("check here")

        auto clipper = TrimmedSurfaceClipping(rAccuracy, rUnit);

        clipper.Clear();

        for (int l = 0; l < mTrimmingLoops.size(); ++l)
        {
            clipper.BeginLoop();

            auto trimming_curves = mTrimmingLoops[l].GetTrimmingCurves();
            //check directions of trimming curves
            if (trimming_curves.size() > 1)
            {
                array_1d<double, 2> end_point_0 = trimming_curves[0].GetCurve2D()->PointAt(
                    trimming_curves[0].GetCurve2D()->Domain().T1());
                array_1d<double, 2> start_point_1 = trimming_curves[1].GetCurve2D()->PointAt(
                    trimming_curves[1].GetCurve2D()->Domain().T0());

                double distance_0_1 = norm_2(end_point_0 - start_point_1);

                array_1d<double, 2> start_point_last = trimming_curves[1].GetCurve2D()->PointAt(
                    trimming_curves[trimming_curves.size() - 1].GetCurve2D()->Domain().T0());

                double distance_0_end = norm_2(end_point_0 - start_point_last);

                if (distance_0_1 < rAccuracy)
                {
                    for (int t = 0; t < trimming_curves.size(); ++t)
                    {
                        auto crv = (trimming_curves[t].GetCurve2D());

                        clipper.AddCurve(*crv);
                    }
                }
                else if (distance_0_end < rAccuracy)
                {
                    for (int t = trimming_curves.size() - 1; t >= 0; --t)
                    {
                        auto crv = (trimming_curves[t].GetCurve2D());

                        clipper.AddCurve(*crv);
                    }
                }
                else
                {
                    KRATOS_ERROR << "Trimming loops of Brep Face " << Id() << " are not defined correctly!" << std::endl;
                }
            }
            else
            {
                for (int t = 0; t < trimming_curves.size(); ++t)
                {
                    auto crv = (trimming_curves[t].GetCurve2D());

                    clipper.AddCurve(*crv);
                }
            }
            clipper.EndLoop();
        }

        clipper.Compute(mNodeSurfaceGeometry3D->SpansU(), mNodeSurfaceGeometry3D->SpansV());

        return clipper;
    }

    //void BrepFace::GetGeometryIntegrationTrimmed(
    //    ModelPart& rModelPart,
    //    const std::string& rType,
    //    const std::string& rName,
    //    const int rShapeFunctionDerivativesOrder,
    //    std::vector<std::string> rVariables) const
    //{
    //    auto clipper = TrimmedSurfaceClipping(0.01, 0.0001);

    //    clipper.Clear();

    //    for (int l = 0; l < mTrimmingLoops.size(); ++l)
    //    {
    //        clipper.BeginLoop();

    //        auto trimming_curves = mTrimmingLoops[l].GetTrimmingCurves();
    //        for (int t = 0; t < trimming_curves.size(); ++t)
    //        {
    //            auto crv = (trimming_curves[t].GetCurve2D());
    //            clipper.AddCurve(*crv);
    //        }

    //        clipper.EndLoop();
    //    }

    //    clipper.Compute(mNodeSurfaceGeometry3D->SpansU(), mNodeSurfaceGeometry3D->SpansV());

    //    int degree_u = mNodeSurfaceGeometry3D->DegreeU();
    //    int degree_v = mNodeSurfaceGeometry3D->DegreeV();

    //    int degree = std::max(degree_u, degree_v) + 1;

    //    for (int i = 0; i < clipper.NbSpansU(); ++i)
    //    {
    //        for (int j = 0; j < clipper.NbSpansU(); ++j)
    //        {
    //            if (clipper.SpanTrimType(i, j) == ANurbs::Empty)
    //            {
    //                continue;
    //            }
    //            if (clipper.SpanTrimType(i, j) == ANurbs::Full)
    //            {
    //                auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
    //                    degree_u + 1,
    //                    degree_v + 1,
    //                    clipper.SpanU(i),
    //                    clipper.SpanV(j));

    //                CreateIntegrationElementsConditions(
    //                    integration_points, rModelPart, rType, rName,
    //                    rShapeFunctionDerivativesOrder, rVariables);

    //                continue;
    //            }
    //            if (clipper.SpanTrimType(i, j) == ANurbs::Trimmed)
    //            {
    //                auto polygons = clipper.SpanPolygons(i, j);

    //                for (int p = 0; p < polygons.size(); ++p)
    //                {
    //                    auto integration_point_polygon = ANurbs::PolygonIntegrationPoints<Kratos::array_1d<double, 2>>();

    //                    integration_point_polygon.Compute(degree, polygons[p]);
    //                    std::vector<ANurbs::IntegrationPoint2<double>> integration_points(integration_point_polygon.NbIntegrationPoints());
    //                    for (int integration_point = 0; integration_point < integration_point_polygon.NbIntegrationPoints(); ++integration_point)
    //                    {
    //                        integration_points[integration_point] = integration_point_polygon.IntegrationPoint(integration_point);
    //                    }
    //                    CreateIntegrationElementsConditions(integration_points, rModelPart, rType, rName,
    //                        rShapeFunctionDerivativesOrder, rVariables);
    //                }
    //            }
    //        }
    //    }
    //}

    void BrepFace::GetIntegrationBrepEdge(
        ModelPart& rModelPart,
        const int trim_index,
        const std::string& rType,
        const std::string& rName,
        const int rShapeFunctionDerivativesOrder,
        std::vector<std::string> rVariables)
    {
        int number_of_points_spans =
            mNodeSurfaceGeometry3D->DegreeU()
            + mNodeSurfaceGeometry3D->DegreeV() + 1;

        auto spans_u = mNodeSurfaceGeometry3D->SpansU();
        auto spans_v = mNodeSurfaceGeometry3D->SpansV();

        ANurbs::SurfaceShapeEvaluator<double> shape(
            mNodeSurfaceGeometry3D->DegreeU(),
            mNodeSurfaceGeometry3D->DegreeV(),
            rShapeFunctionDerivativesOrder);

        auto curve_2d = GetTrimCurve(trim_index);

        auto curve_on_surface_3d = Kratos::make_shared<CurveOnSurface<3>>(
                curve_2d->CurveGeometry(), mNodeSurfaceGeometry3D, curve_2d->Domain());

        auto curve_knot_intersections = curve_on_surface_3d->Spans();

        for (int i = 0; i < curve_knot_intersections.size(); ++i)
        {
            auto integration_points = ANurbs::IntegrationPoints<double>::Points1(
                number_of_points_spans, curve_knot_intersections[i]);

            for (int ip = 0; ip < integration_points.size(); ++ip)
            {
                auto derivatives = curve_2d->DerivativesAt(integration_points[ip].t, 1);

                Element::GeometryType::PointsArrayType control_points;
                Vector shape_function;
                Matrix shape_function_derivative;
                Matrix shape_function_second_derivative;

                EvaluatePoint(derivatives[0][0], derivatives[0][1],
                    control_points,
                    shape_function,
                    shape_function_derivative,
                    shape_function_second_derivative);

                if (rType == "element")
                {
                    int id = 0;
                    if (rModelPart.GetRootModelPart().Elements().size() > 0)
                        id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

                    auto element = rModelPart.CreateNewElement(rName, id, control_points, 0);

                    rModelPart.AddNodes(control_points.begin(), control_points.end());


                    element->SetValue(SHAPE_FUNCTION_VALUES, shape_function);
                    element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, shape_function_derivative);
                    element->SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, shape_function_second_derivative);

                    Vector tangents(2);
                    tangents[0] = derivatives[1][0];
                    tangents[1] = derivatives[1][1];
                    element->SetValue(TANGENTS, tangents);

                    element->SetValue(INTEGRATION_WEIGHT, integration_points[ip].weight);
                    element->SetValue(BREP_ID, this->Id());
                }
                if (rType == "condition")
                {
                    int id = 0;
                    if (rModelPart.GetRootModelPart().Conditions().size() > 0)
                        id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;

                    auto condition = rModelPart.CreateNewCondition(rName, id, control_points, 0);

                    rModelPart.AddNodes(control_points.begin(), control_points.end());


                    condition->SetValue(SHAPE_FUNCTION_VALUES, shape_function);
                    condition->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, shape_function_derivative);
                    condition->SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, shape_function_second_derivative);

                    Vector tangents(2);
                    tangents[0] = derivatives[1][0];
                    tangents[1] = derivatives[1][1];
                    condition->SetValue(TANGENTS, tangents);

                    condition->SetValue(INTEGRATION_WEIGHT, integration_points[ip].weight);
                    condition->SetValue(BREP_ID, this->Id());
                }
                // for strong application of properties on control point nodes
                if (rType == "node")
                {
                    for (int sh_nonzero = 0; sh_nonzero < shape_function.size(); ++sh_nonzero)
                    {
                        if (shape_function(sh_nonzero) > 1e-5)
                        {
                            rModelPart.AddNode(control_points(sh_nonzero));
                        }
                    }
                }
            }
        }
    }

    //void BrepFace::CreateIntegrationElementsConditions(
    //    std::vector<ANurbs::IntegrationPoint2<double>> rIntegrationPoints,
    //    ModelPart& rModelPart,
    //    const std::string& rType,
    //    const std::string& rName,
    //    const int rShapeFunctionDerivativesOrder,
    //    std::vector<std::string> rVariables) const
    //{
    //    ANurbs::SurfaceShapeEvaluator<double> shape(
    //        mNodeSurfaceGeometry3D->DegreeU(),
    //        mNodeSurfaceGeometry3D->DegreeV(),
    //        rShapeFunctionDerivativesOrder);

    //    for (int k = 0; k < rIntegrationPoints.size(); ++k)
    //    {
    //        array_1d<double, 2> local_coordinates;
    //        local_coordinates[0] = rIntegrationPoints[k].u;
    //        local_coordinates[1] = rIntegrationPoints[k].v;

    //        shape.Compute(
    //            mNodeSurfaceGeometry3D->KnotsU(),
    //            mNodeSurfaceGeometry3D->KnotsV(),
    //            mNodeSurfaceGeometry3D->Weights(),
    //            rIntegrationPoints[k].u,
    //            rIntegrationPoints[k].v);

    //        Element::GeometryType::PointsArrayType non_zero_control_points;

    //        Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
    //        Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 2);
    //        Matrix N_2 = ZeroMatrix(shape.NbNonzeroPoles(), 3);

    //        for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
    //        {
    //            int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
    //            int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

    //            non_zero_control_points.push_back(mNodeSurfaceGeometry3D->GetNode(
    //                shape.NonzeroPoleIndices()[n].first, shape.NonzeroPoleIndices()[n].second));

    //            N_0[n] = shape(0, indexU, indexV);
    //            N_1(n, 0) = shape(1, indexU, indexV);
    //            N_1(n, 1) = shape(2, indexU, indexV);
    //            N_2(n, 0) = shape(3, indexU, indexV);
    //            N_2(n, 1) = shape(5, indexU, indexV);
    //            N_2(n, 2) = shape(4, indexU, indexV);
    //        }


    //        if (rType == "element")
    //        {
    //            int id = 0;
    //            if (rModelPart.GetRootModelPart().Elements().size() > 0)
    //                id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

    //            rModelPart.AddNodes(non_zero_control_points.begin(), non_zero_control_points.end());

    //            auto element = rModelPart.CreateNewElement(rName, id, non_zero_control_points, 0);

    //            element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
    //            element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
    //            element->SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, N_2);
    //            element->SetValue(INTEGRATION_WEIGHT, rIntegrationPoints[k].weight);

    //            element->SetValue(LOCAL_COORDINATES, local_coordinates);
    //            element->SetValue(BREP_ID, this->Id());
    //        }

    //        if (rType == "condition")
    //        {
    //            int id = 0;
    //            if (rModelPart.GetRootModelPart().Conditions().size() > 0)
    //                id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;

    //            rModelPart.AddNodes(non_zero_control_points.begin(), non_zero_control_points.end());

    //            auto condition = rModelPart.CreateNewCondition(rName, id, non_zero_control_points, 0);

    //            condition->SetValue(SHAPE_FUNCTION_VALUES, N_0);
    //            condition->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
    //            condition->SetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, N_2);
    //            condition->SetValue(INTEGRATION_WEIGHT, rIntegrationPoints[k].weight);

    //            condition->SetValue(LOCAL_COORDINATES, local_coordinates);
    //            condition->SetValue(BREP_ID, this->Id());
    //        }

    //        // for strong application of properties on control point nodes
    //        if (rType == "node")
    //        {
    //            for (int sh_nonzero = 0; sh_nonzero < N_0.size(); ++sh_nonzero)
    //            {
    //                if (N_0(sh_nonzero) > 1e-5)
    //                {
    //                    rModelPart.AddNode(non_zero_control_points(sh_nonzero));
    //                }
    //            }
    //        }
    //    }
    //}

    void BrepFace::EvaluatePoint(
        const double rU,
        const double rV,
        Element::GeometryType::PointsArrayType& rNonZeroControlPoints,
        Vector& rShapeFunction,
        Matrix& rShapeFunctionDerivative,
        Matrix& rShapeFunctionSecondDerivative
    ) const
    {
        ANurbs::SurfaceShapeEvaluator<double> shape(
            mNodeSurfaceGeometry3D->DegreeU(),
            mNodeSurfaceGeometry3D->DegreeV(),
            3);

        shape.Compute(
            mNodeSurfaceGeometry3D->KnotsU(),
            mNodeSurfaceGeometry3D->KnotsV(),
            mNodeSurfaceGeometry3D->Weights(),
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
        array_1d<double,3> location3D = ZeroVector(3);
        for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
        {
            int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
            int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

            rNonZeroControlPoints.push_back(mNodeSurfaceGeometry3D->GetNode(
                shape.NonzeroPoleIndices()[n].first, shape.NonzeroPoleIndices()[n].second));

            rShapeFunction[n] = shape(0, indexU, indexV);
            rShapeFunctionDerivative(n, 0) = shape(1, indexU, indexV);
            rShapeFunctionDerivative(n, 1) = shape(2, indexU, indexV);
            rShapeFunctionSecondDerivative(n, 0) = shape(3, indexU, indexV);
            rShapeFunctionSecondDerivative(n, 1) = shape(5, indexU, indexV);
            rShapeFunctionSecondDerivative(n, 2) = shape(4, indexU, indexV);
        }
    }

    const Kratos::shared_ptr<Curve<2>> BrepFace::GetTrimCurve(const int trim_index) const
    {
        for (int i = 0; i < mTrimmingLoops.size(); ++i)
        {
            std::vector<BrepTrimmingCurve> trimming_curves = mTrimmingLoops[i].GetTrimmingCurves();
            for (int j = 0; j < trimming_curves.size(); ++j)
            {
                if (trimming_curves[j].GetTrimIndex() == trim_index)
                    return trimming_curves[j].GetCurve2D();
            }
        }
        for (int i = 0; i < mEmbeddedLoops.size(); ++i)
        {
            std::vector<BrepTrimmingCurve> trimming_loops = mTrimmingLoops[i].GetTrimmingCurves();
            for (int j = 0; j < trimming_loops.size(); ++j)
            {
                if (trimming_loops[j].GetTrimIndex() == trim_index)
                    return trimming_loops[j].GetCurve2D();
            }
        }
        for (int i = 0; i < mEmbeddedEdges.size(); ++i)
        {
            if (mEmbeddedEdges[i].GetTrimIndex() == trim_index)
                return mEmbeddedEdges[i].GetCurve2D();
        }
        KRATOS_ERROR << "Trimming curve of index: " << trim_index << " does not exist in BrepFace " << Id() << std::endl;
    }

    const Kratos::shared_ptr<NodeSurfaceGeometry3D> BrepFace::GetSurface() const
    {
        return mNodeSurfaceGeometry3D;
    }
    

    const bool BrepFace::GetIsTrimmed() const 
    {
        return m_is_trimmed; 
    }

    const bool BrepFace::GetIsRational() const
    {
        return m_is_rational; 
    }

    const std::vector<BrepBoundaryLoop> BrepFace::GetBoundaryLoop() const
    { 
        return mTrimmingLoops;      
    }


    ///Constructor
    BrepFace::BrepFace(
        int rBrepId,
        bool rIsTrimmed,
        bool rIsRational,
        std::vector<BrepBoundaryLoop>& rTrimmingLoops,
        std::vector<BrepBoundaryLoop>& rEmbeddedLoops,
        std::vector<BrepTrimmingCurve>& rEmbeddedEdges,
        std::vector<EmbeddedPoint>& rEmbeddedPoints,
        Vector& rKnotVectorU,
        Vector& rKnotVectorV,
        int rP,
        int rQ,
        IntVector& rControlPointIds,
        ModelPart& rModelPart)
        : mTrimmingLoops(rTrimmingLoops),
          m_is_trimmed(rIsTrimmed),
          m_is_rational(rIsRational),
          mEmbeddedLoops(rEmbeddedLoops),
          mEmbeddedEdges(rEmbeddedEdges),
          mEmbeddedPoints(rEmbeddedPoints),
          mModelPart(rModelPart),
          IndexedObject(rBrepId),
          Flags()
    {
        int number_of_nodes_u = rKnotVectorU.size() - rP - 1;
        int number_of_nodes_v = rKnotVectorV.size() - rQ - 1;
        
        mNodeSurfaceGeometry3D = Kratos::make_unique<NodeSurfaceGeometry3D>(
            rP, rQ, number_of_nodes_u, number_of_nodes_v);

        for (int i = 0; i < rKnotVectorU.size()-2; ++i)
        {
            mNodeSurfaceGeometry3D->SetKnotU(
                i, rKnotVectorU(i+1));
        }

        for (int i = 0; i < rKnotVectorV.size()-2; ++i)
        {
            mNodeSurfaceGeometry3D->SetKnotV(i, rKnotVectorV(i+1));
        }

        for (int i = 0; i < number_of_nodes_u; ++i)
        {
            for (int j = 0; j < number_of_nodes_v; ++j)
            {
                Node<3>::Pointer node = rModelPart.pGetNode(
                    rControlPointIds[j*(number_of_nodes_u)+i]);
                mNodeSurfaceGeometry3D->SetNode(i, j, node);
                if (rIsRational)
                {
                    mNodeSurfaceGeometry3D->SetWeight(i, j,
                        node->GetValue(NURBS_CONTROL_POINT_WEIGHT));
                }
            }
        }
    }
} // namespace Kratos.