//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

#if !defined(IGA_INTEGRATION_UTILITIES_H_INCLUDED)
#define IGA_INTEGRATION_UTILITIES_H_INCLUDED

#include "includes/model_part.h"
#include "geometries/integration_point_curve_on_surface_3d.h"
#include "geometries/integration_point_surface_3d.h"
#include "geometries/integration_point_point_on_surface_3d.h"
#include "geometries/coupling_master_slave.h"

namespace Kratos
{

namespace IgaIntegrationUtilities
{
    static void ChangeElementType(
        std::vector<Element::Pointer> rElementList,
        ModelPart& rModelPart,
        std::string& rElementName,
        int& rIdCounter)
    {
        const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

        ModelPart::ElementsContainerType new_element_list;
        new_element_list.reserve(rElementList.size());

        for (std::size_t i = 0; i< rElementList.size(); i++) 
        {
            auto it_elem = rElementList[i];

            auto p_element = rReferenceElement.Create(rIdCounter, it_elem->pGetGeometry(), nullptr);

            rIdCounter++;

            // Deep copy elemental data and flags
            p_element->Data() = it_elem->Data();
            p_element->Set(Flags(*it_elem));

            rModelPart.AddNodes(p_element->GetGeometry().begin(), p_element->GetGeometry().end());

            new_element_list.push_back(p_element);
        }

        rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }
    static void ChangeConditionType(
        std::vector<Element::Pointer> rElementList,
        ModelPart& rModelPart,
        std::string& rConditionName,
        int& rIdCounter)
    {
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;
        new_condition_list.reserve(rElementList.size());

        for (std::size_t i = 0; i< rElementList.size(); i++)
        {
            auto it_elem = rElementList[i];

            auto p_condition = rReferenceCondition.Create(rIdCounter, it_elem->pGetGeometry(), nullptr);

            rIdCounter++;

            // Deep copy elemental data and flags
            p_condition->Data() = it_elem->Data();
            p_condition->Set(Flags(*it_elem));

            new_condition_list.push_back(p_condition);
        }

        rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
    }

    static void CreateElements(
        std::vector<Geometry<Node<3>>::Pointer> rGeometryList,
        ModelPart& rModelPart,
        std::string& rElementName,
        int& rIdCounter)
    {
        const Element& rReferenceElement = KratosComponents<Element>::Get(rElementName);

        ModelPart::ElementsContainerType new_element_list;
        new_element_list.reserve(rGeometryList.size());

        for (int i = 0; i < static_cast<int>(rGeometryList.size()); i++)
        {
            auto p_element = rReferenceElement.Create(rIdCounter, rGeometryList[i], nullptr);

            rIdCounter++;

            rModelPart.AddNodes(rGeometryList[i]->begin(), rGeometryList[i]->end());

            new_element_list.push_back(p_element);
        }

        rModelPart.AddElements(new_element_list.begin(), new_element_list.end());
    }

    static void CreateConditions(
        std::vector<Geometry<Node<3>>::Pointer> rGeometryList,
        ModelPart& rModelPart,
        std::string& rConditionName,
        int& rIdCounter)
    {
        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;
        new_condition_list.reserve(rGeometryList.size());

        for (int i = 0; i< static_cast<int>(rGeometryList.size()); i++)
        {
            auto p_condition = rReferenceCondition.Create(rIdCounter, rGeometryList[i], nullptr);

            rIdCounter++;

            rModelPart.AddNodes(rGeometryList[i]->begin(), rGeometryList[i]->end());

            new_condition_list.push_back(p_condition);
        }

        rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());
    }


    static std::vector<Element::Pointer> GetIntegrationDomainGeometrySurface(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
        const TrimmedSurfaceClipping& rClipper,
        int ShapeFunctionDerivativesOrder)
    {
        std::vector<Element::Pointer> new_elements;

        int degree_u = pSurface->DegreeU();
        int degree_v = pSurface->DegreeV();

        int degree = std::max(degree_u, degree_v) + 1;

        ANurbs::SurfaceShapeEvaluator<double> shape(
            pSurface->DegreeU(),
            pSurface->DegreeV(),
            ShapeFunctionDerivativesOrder);

        for (int i = 0; i < rClipper.NbSpansU(); ++i)
        {
            for (int j = 0; j < rClipper.NbSpansV(); ++j)
            {
                if (rClipper.SpanTrimType(i, j) == ANurbs::Empty)
                {
                    continue;
                }
                else if (rClipper.SpanTrimType(i, j) == ANurbs::Full)
                {
                    auto integration_points = ANurbs::IntegrationPoints<double>::Points2(
                        degree_u + 1,
                        degree_v + 1,
                        rClipper.SpanU(i),
                        rClipper.SpanV(j));

                    for (int i = 0; i < integration_points.size(); ++i)
                    {
                        array_1d<double, 2> local_coordinates;
                        local_coordinates[0] = integration_points[i].u;
                        local_coordinates[1] = integration_points[i].v;

                        shape.Compute(
                            pSurface->KnotsU(),
                            pSurface->KnotsV(),
                            pSurface->Weights(),
                            integration_points[i].u,
                            integration_points[i].v);

                        Element::GeometryType::PointsArrayType non_zero_control_points;
                        Vector shape_function(shape.NonzeroPoleIndices().size());
                        Matrix shape_function_derivative(shape.NonzeroPoleIndices().size(), 2);
                        Matrix shape_function_second_derivative(shape.NonzeroPoleIndices().size(), 3);

                        for (int n = 0; n < shape.NonzeroPoleIndices().size(); ++n)
                        {
                            int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                            int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                            non_zero_control_points.push_back(pSurface->GetNode(
                                shape.NonzeroPoleIndices()[n].first,
                                shape.NonzeroPoleIndices()[n].second));

                            if (ShapeFunctionDerivativesOrder > -1)
                                shape_function[n] = shape(0, indexU, indexV);
                            if (ShapeFunctionDerivativesOrder > 0)
                            {
                                shape_function_derivative(n, 0) = shape(1, indexU, indexV);
                                shape_function_derivative(n, 1) = shape(2, indexU, indexV);
                            }
                            if (ShapeFunctionDerivativesOrder > 1)
                            {
                                shape_function_second_derivative(n, 0) = shape(3, indexU, indexV);
                                shape_function_second_derivative(n, 1) = shape(5, indexU, indexV);
                                shape_function_second_derivative(n, 2) = shape(4, indexU, indexV);
                            }
                        }

                        Element ele(0, non_zero_control_points);

                        Element::Pointer element = Kratos::make_intrusive<Element>(ele);

                        if (ShapeFunctionDerivativesOrder > -1)
                        {
                            element->SetValue(
                                SHAPE_FUNCTION_VALUES,
                                shape_function);
                        }
                        if (ShapeFunctionDerivativesOrder > 0)
                        {
                            element->SetValue(
                                SHAPE_FUNCTION_LOCAL_DERIVATIVES,
                                shape_function_derivative);
                        }
                        if (ShapeFunctionDerivativesOrder > 1)
                        {
                            element->SetValue(
                                SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES,
                                shape_function_second_derivative);
                        }
                        element->SetValue(
                            INTEGRATION_WEIGHT,
                            integration_points[i].weight);

                        element->SetValue(LOCAL_COORDINATES, local_coordinates);

                        new_elements.push_back(element);
                    }
                }
                else if (rClipper.SpanTrimType(i, j) == ANurbs::Trimmed)
                {
                    auto polygons = rClipper.SpanPolygons(i, j);

                    for (int p = 0; p < polygons.size(); ++p)
                    {
                        auto integration_point_polygon = ANurbs::PolygonIntegrationPoints<Kratos::array_1d<double, 2>>();

                        integration_point_polygon.Compute(degree, polygons[p]);

                        for (int i = 0; i < integration_point_polygon.NbIntegrationPoints(); ++i)
                        {
                            array_1d<double, 2> local_coordinates;
                            local_coordinates[0] = integration_point_polygon.IntegrationPoint(i).u;
                            local_coordinates[1] = integration_point_polygon.IntegrationPoint(i).v;

                            shape.Compute(
                                pSurface->KnotsU(),
                                pSurface->KnotsV(),
                                pSurface->Weights(),
                                integration_point_polygon.IntegrationPoint(i).u,
                                integration_point_polygon.IntegrationPoint(i).v);

                            int number_of_non_zero_cps = shape.NonzeroPoleIndices().size();
                            Element::GeometryType::PointsArrayType non_zero_control_points;
                            Vector shape_function(number_of_non_zero_cps);
                            Matrix shape_function_derivative(number_of_non_zero_cps, 2);
                            Matrix shape_function_second_derivative(number_of_non_zero_cps, 3);

                            for (int n = 0; n < number_of_non_zero_cps; ++n)
                            {
                                int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                                int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                                non_zero_control_points.push_back(pSurface->GetNode(
                                    shape.NonzeroPoleIndices()[n].first,
                                    shape.NonzeroPoleIndices()[n].second));

                                if (ShapeFunctionDerivativesOrder > -1)
                                    shape_function[n] = shape(0, indexU, indexV);
                                if (ShapeFunctionDerivativesOrder > 0)
                                {
                                    shape_function_derivative(n, 0) = shape(1, indexU, indexV);
                                    shape_function_derivative(n, 1) = shape(2, indexU, indexV);
                                }
                                if (ShapeFunctionDerivativesOrder > 1)
                                {
                                    shape_function_second_derivative(n, 0) = shape(3, indexU, indexV);
                                    shape_function_second_derivative(n, 1) = shape(5, indexU, indexV);
                                    shape_function_second_derivative(n, 2) = shape(4, indexU, indexV);
                                }
                            }

                            Element ele(0, non_zero_control_points);

                            Element::Pointer element = Kratos::make_intrusive<Element>(ele);

                            if (ShapeFunctionDerivativesOrder > -1)
                            {
                                element->SetValue(
                                    SHAPE_FUNCTION_VALUES,
                                    shape_function);
                            }
                            if (ShapeFunctionDerivativesOrder > 0)
                            {
                                element->SetValue(
                                    SHAPE_FUNCTION_LOCAL_DERIVATIVES,
                                    shape_function_derivative);
                            }
                            if (ShapeFunctionDerivativesOrder > 1)
                            {
                                element->SetValue(
                                    SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES,
                                    shape_function_second_derivative);
                            }
                            element->SetValue(
                                INTEGRATION_WEIGHT,
                                integration_point_polygon.IntegrationPoint(i).weight);

                            element->SetValue(LOCAL_COORDINATES, local_coordinates);

                            new_elements.push_back(element);
                        }
                    }
                }
            }
        }
        return new_elements;
    }

    static std::vector<Geometry<Node<3>>::Pointer> GetIntegrationDomainSurfaceEdge(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
        const std::shared_ptr<Curve<2>>& pTrimmingCurve,
        int ShapeFunctionDerivativesOrder)
    {
        std::vector<Geometry<Node<3>>::Pointer> new_geometries;

        int number_of_points_spans =
            pSurface->DegreeU() + pSurface->DegreeV() + 1;

        auto spans_u = pSurface->SpansU();
        auto spans_v = pSurface->SpansV();

        ANurbs::SurfaceShapeEvaluator<double> shape(
            pSurface->DegreeU(),
            pSurface->DegreeV(),
            ShapeFunctionDerivativesOrder);

        auto curve_on_surface_3d = std::make_shared<CurveOnSurface<3>>(
            pTrimmingCurve->CurveGeometry(), pSurface, pTrimmingCurve->Domain());

        auto curve_knot_intersections = curve_on_surface_3d->Spans();

        for (int i = 0; i < curve_knot_intersections.size(); ++i)
        {
            auto integration_points = ANurbs::IntegrationPoints<double>::Points1(
                number_of_points_spans, curve_knot_intersections[i]);

            for (int ip = 0; ip < integration_points.size(); ++ip)
            {
                auto derivatives = pTrimmingCurve->DerivativesAt(integration_points[ip].t, 1);

                shape.Compute(
                    pSurface->KnotsU(),
                    pSurface->KnotsV(),
                    pSurface->Weights(),
                    derivatives[0][0],
                    derivatives[0][1]);

                Element::GeometryType::PointsArrayType control_points;
                Vector shape_function;
                Matrix shape_function_derivative;
                Matrix shape_function_second_derivative;

                Geometry<Node<3>>::PointsArrayType cps;

                int number_of_non_zero_cps = shape.NonzeroPoleIndices().size();

                Matrix N(1, number_of_non_zero_cps);
                Matrix DN_De(number_of_non_zero_cps, 2);

                for (int n = 0; n < number_of_non_zero_cps; ++n)
                {
                    int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                    int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                    cps.push_back(pSurface->GetNode(
                        shape.NonzeroPoleIndices()[n].first,
                        shape.NonzeroPoleIndices()[n].second));

                    if (ShapeFunctionDerivativesOrder > -1)
                        N(0, n) = shape(0, indexU, indexV);

                    if (ShapeFunctionDerivativesOrder > 0)
                    {
                        DN_De(n, 0) = shape(1, indexU, indexV);
                        DN_De(n, 1) = shape(2, indexU, indexV);
                    }
                }

                array_1d<double, 2> local_coordinates;
                local_coordinates[0] = derivatives[0][0];
                local_coordinates[1] = derivatives[0][1];

                array_1d<double, 2> tangents;
                tangents[0] = derivatives[1][0];
                tangents[1] = derivatives[1][1];

                Geometry<Node<3>>::ShapeFunctionsGradientsType DN_De_type(1);
                DN_De_type[0] = DN_De;

                Geometry<Node<3>>::IntegrationPointsArrayType ips(1);
                ips[0] = IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointType(local_coordinates[0], local_coordinates[1], integration_points[ip].weight);

                IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointsContainerType ips_container =
                { { ips } };
                IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsValuesContainerType N_container =
                { { N } };
                IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsLocalGradientsContainerType DN_De_container =
                { { DN_De_type } };

                Geometry<Node<3>>::Pointer intgeration_point_geometry =
                    std::make_shared<IntegrationPointCurveOnSurface3d<Node<3>>>(
                        cps,
                        ips_container,
                        tangents,
                        N_container,
                        DN_De_container);

                new_geometries.push_back(intgeration_point_geometry);
            }
        }

        return new_geometries;
    }

    static Geometry<Node<3>>::Pointer GetIntegrationPointSurfacePoint(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
        const double U,
        const double V,
        int ShapeFunctionDerivativesOrder)
    {
        ANurbs::SurfaceShapeEvaluator<double> shape(
            pSurface->DegreeU(),
            pSurface->DegreeV(),
            ShapeFunctionDerivativesOrder);

        shape.Compute(
            pSurface->KnotsU(),
            pSurface->KnotsV(),
            pSurface->Weights(),
            U,
            V);

        Geometry<Node<3>>::PointsArrayType cps;

        int number_of_non_zero_cps_1 = shape.NonzeroPoleIndices().size();

        Matrix shape_function(1, number_of_non_zero_cps_1);
        Matrix shape_function_derivative(number_of_non_zero_cps_1, 2);
        //Matrix shape_function_second_derivative(number_of_non_zero_cps_1,3);
        //KRATOS_WATCH(number_of_non_zero_cps_1)

        for (int n = 0; n < number_of_non_zero_cps_1; ++n)
        {
            int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
            int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

            cps.push_back(pSurface->GetNode(
                shape.NonzeroPoleIndices()[n].first,
                shape.NonzeroPoleIndices()[n].second));

            if (ShapeFunctionDerivativesOrder > -1)
                shape_function(0, n) = shape(0, indexU, indexV);
            if (ShapeFunctionDerivativesOrder > 0)
            {
                shape_function_derivative(n, 0) = shape(1, indexU, indexV);
                shape_function_derivative(n, 1) = shape(2, indexU, indexV);
            }
            //if (ShapeFunctionDerivativesOrder > 1)
            //{
            //    shape_function_second_derivative(n, 0) = shape_1(3, indexU, indexV);
            //    shape_function_second_derivative(n, 1) = shape_1(5, indexU, indexV);
            //    shape_function_second_derivative(n, 2) = shape_1(4, indexU, indexV);
            //}
        }

        Geometry<Node<3>>::ShapeFunctionsGradientsType shape_function_gradients(1);
        shape_function_gradients[0] = shape_function_derivative;

        Geometry<Node<3>>::IntegrationPointsArrayType ips(1);
        ips[0] = IntegrationPointPointOnSurface3d<Node<3>>::IntegrationPointType(
            U,
            V,
            1.0);

        IntegrationPointPointOnSurface3d<Node<3>>::IntegrationPointsContainerType IntegrationPointsArray =
        { { ips } };
        IntegrationPointPointOnSurface3d<Node<3>>::ShapeFunctionsValuesContainerType sfcontainer =
        { { shape_function } };
        IntegrationPointPointOnSurface3d<Node<3>>::ShapeFunctionsLocalGradientsContainerType sfgcontainer =
        { { shape_function_gradients } };


        Geometry<Node<3>>::Pointer intgeration_point_point =
            Kratos::make_shared<IntegrationPointPointOnSurface3d<Node<3>>>(
            cps,
            IntegrationPointsArray,
            sfcontainer,
            sfgcontainer);

        return intgeration_point_point;
    }

    static std::vector<Geometry<Node<3>>::Pointer> GetIntegrationDomainSurfaceEdgeSurfaceEdge(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface1,
        const Kratos::shared_ptr<Curve<2>>& pTrimmingCurve1,
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface2,
        const Kratos::shared_ptr<Curve<2>>& pTrimmingCurve2,
        int ShapeFunctionDerivativesOrder)
    {
        std::vector<Geometry<Node<3>>::Pointer> new_geometries;

        int degree_1 = pSurface1->DegreeU() + pSurface1->DegreeV();
        int degree_2 = pSurface2->DegreeU() + pSurface2->DegreeV();

        int number_of_points_per_knot_span = std::max(degree_1, degree_2) + 1;

        auto curve_on_surface_3d_1 = Kratos::make_shared<CurveOnSurface<3>>(
            pTrimmingCurve1->CurveGeometry(),
            pSurface1,
            pTrimmingCurve1->Domain());
        auto curve_on_surface_3d_2 = Kratos::make_shared<CurveOnSurface<3>>(
            pTrimmingCurve2->CurveGeometry(),
            pSurface2,
            pTrimmingCurve2->Domain());

        ANurbs::SurfaceShapeEvaluator<double> shape_1(
            pSurface1->DegreeU(),
            pSurface1->DegreeV(),
            ShapeFunctionDerivativesOrder);

        ANurbs::SurfaceShapeEvaluator<double> shape_2(
            pSurface2->DegreeU(),
            pSurface2->DegreeV(),
            ShapeFunctionDerivativesOrder);

        auto projection_1 = ANurbs::PointOnCurveProjection<Kratos::array_1d<double, 3>>(
            curve_on_surface_3d_1,
            0.1);
        auto projection_2 = ANurbs::PointOnCurveProjection<Kratos::array_1d<double, 3>>(
            curve_on_surface_3d_2,
            0.1);

        auto curve_knot_intersections_2 = curve_on_surface_3d_2->Spans();

        std::vector<double> curve_knot_intersections_vector;
        for (int i = 0; i < curve_knot_intersections_2.size(); ++i)
        {
            curve_knot_intersections_vector.push_back(curve_knot_intersections_2[i].T0());
            auto point_3d = curve_on_surface_3d_1->PointAt(curve_knot_intersections_2[i].T0());
            projection_1.Compute(point_3d);
            curve_knot_intersections_vector.push_back(projection_1.Parameter());
        }
        curve_knot_intersections_vector.push_back(curve_on_surface_3d_2->Domain().T1());
        std::sort(curve_knot_intersections_vector.begin(), curve_knot_intersections_vector.end());

        for (int i = 0; i < curve_knot_intersections_vector.size() - 1; ++i)
        {
            if (std::abs(curve_knot_intersections_vector[i] - curve_knot_intersections_vector[i + 1]) > 1e-6)
            {
                ANurbs::Interval<double> interval(curve_knot_intersections_vector[i], curve_knot_intersections_vector[i + 1]);

                auto integration_points = ANurbs::IntegrationPoints<double>::Points1(number_of_points_per_knot_span, interval);

                for (int ip = 0; ip < integration_points.size(); ++ip)
                {
                    //MASTER

                    auto point_3d = curve_on_surface_3d_1->PointAt(integration_points[ip].t);
                    auto derivatives_1 = pTrimmingCurve1->DerivativesAt(integration_points[ip].t, 1);

                    shape_1.Compute(
                        pSurface1->KnotsU(),
                        pSurface1->KnotsV(),
                        pSurface1->Weights(),
                        derivatives_1[0][0],
                        derivatives_1[0][1]);

                    Geometry<Node<3>>::PointsArrayType cps_master;

                    int number_of_non_zero_cps_1 = shape_1.NonzeroPoleIndices().size();

                    Matrix shape_function_master(1, number_of_non_zero_cps_1);
                    Matrix shape_function_derivative_master(number_of_non_zero_cps_1, 2);
                    //Matrix shape_function_second_derivative(number_of_non_zero_cps_1,3);
                    //KRATOS_WATCH(number_of_non_zero_cps_1)

                    for (int n = 0; n < number_of_non_zero_cps_1; ++n)
                    {
                        int indexU = shape_1.NonzeroPoleIndices()[n].first - shape_1.FirstNonzeroPoleU();
                        int indexV = shape_1.NonzeroPoleIndices()[n].second - shape_1.FirstNonzeroPoleV();

                        cps_master.push_back(pSurface1->GetNode(
                            shape_1.NonzeroPoleIndices()[n].first,
                            shape_1.NonzeroPoleIndices()[n].second));

                        if (ShapeFunctionDerivativesOrder>-1)
                            shape_function_master(0, n) = shape_1(0, indexU, indexV);
                        if (ShapeFunctionDerivativesOrder > 0)
                        {
                            shape_function_derivative_master(n, 0) = shape_1(1, indexU, indexV);
                            shape_function_derivative_master(n, 1) = shape_1(2, indexU, indexV);
                        }
                        //if (ShapeFunctionDerivativesOrder > 1)
                        //{
                        //    shape_function_second_derivative(n, 0) = shape_1(3, indexU, indexV);
                        //    shape_function_second_derivative(n, 1) = shape_1(5, indexU, indexV);
                        //    shape_function_second_derivative(n, 2) = shape_1(4, indexU, indexV);
                        //}
                    }

                    Geometry<Node<3>>::ShapeFunctionsGradientsType shape_function_gradients;
                    //shape_function_gradients.back(shape_function_derivative);

                    array_1d<double, 2> tangents;
                    tangents[0] = derivatives_1[1][0];
                    tangents[1] = derivatives_1[1][1];

                    array_1d<double, 2> local_coordinates;
                    local_coordinates[0] = derivatives_1[0][0];
                    local_coordinates[1] = derivatives_1[0][1];

                    Geometry<Node<3>>::ShapeFunctionsGradientsType shape_function_gradients_master(1);
                    shape_function_gradients_master[0] = shape_function_derivative_master;

                    Geometry<Node<3>>::IntegrationPointsArrayType ips(1);
                    ips[0] = IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointType(local_coordinates[0], local_coordinates[1], integration_points[ip].weight);

                    IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointsContainerType IntegrationPointsArray =
                    { { ips } };
                    IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsValuesContainerType sfcontainer_master =
                    { { shape_function_master } };
                    IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsLocalGradientsContainerType sfgcontainer_master =
                    { { shape_function_gradients_master } };


                    Geometry<Node<3>>::Pointer intgeration_point_master = Kratos::make_shared<IntegrationPointCurveOnSurface3d<Node<3>>>(
                        cps_master,
                        IntegrationPointsArray,
                        tangents,
                        sfcontainer_master,
                        sfgcontainer_master);
                    //Geometry<Node<3>> intgeration_point_master = IntegrationPointCurveOnSurface3d<Node<3>>(
                    //    cps_master,
                    //    local_coordinates,
                    //    tangents,
                    //    integration_points[ip].weight,
                    //    shape_function,
                    //    shape_function_gradients);

                    //SLAVE
                    projection_2.Compute(point_3d);
                    auto derivatives_2 = pTrimmingCurve2->DerivativesAt(projection_2.Parameter(), 1);

                    shape_2.Compute(
                        pSurface2->KnotsU(),
                        pSurface2->KnotsV(),
                        pSurface2->Weights(),
                        derivatives_2[0][0],
                        derivatives_2[0][1]);

                    Geometry<Node<3>>::PointsArrayType cps_slave;
                    int number_of_non_zero_cps_2 = shape_2.NonzeroPoleIndices().size();

                    Matrix shape_function_slave(1, number_of_non_zero_cps_2);
                    Matrix shape_function_derivative_slave(number_of_non_zero_cps_2, 2);

                    for (int n = 0; n < number_of_non_zero_cps_2; ++n)
                    {
                        int indexU = shape_2.NonzeroPoleIndices()[n].first - shape_2.FirstNonzeroPoleU();
                        int indexV = shape_2.NonzeroPoleIndices()[n].second - shape_2.FirstNonzeroPoleV();

                        cps_slave.push_back(pSurface2->GetNode(
                            shape_2.NonzeroPoleIndices()[n].first,
                            shape_2.NonzeroPoleIndices()[n].second));

                        if (ShapeFunctionDerivativesOrder>-1)
                            shape_function_slave(0, n) = shape_2(0, indexU, indexV);
                        if (ShapeFunctionDerivativesOrder > 0)
                        {
                            shape_function_derivative_slave(n, 0) = shape_2(1, indexU, indexV);
                            shape_function_derivative_slave(n, 1) = shape_2(2, indexU, indexV);
                        }
                        //if (ShapeFunctionDerivativesOrder > 1)
                        //{
                        //    shape_function_second_derivative_slave(n, 0) = shape_2(3, indexU, indexV);
                        //    shape_function_second_derivative_slave(n, 1) = shape_2(5, indexU, indexV);
                        //    shape_function_second_derivative_slave(n, 2) = shape_2(4, indexU, indexV);
                        //}
                    }

                    array_1d<double, 2> tangents_slave;
                    tangents_slave[0] = derivatives_2[1][0];
                    tangents_slave[1] = derivatives_2[1][1];

                    array_1d<double, 2> local_coordinates_slave;
                    local_coordinates_slave[0] = derivatives_2[0][0];
                    local_coordinates_slave[1] = derivatives_2[0][1];

                    Geometry<Node<3>>::ShapeFunctionsGradientsType shape_function_gradients_slave(1);
                    shape_function_gradients_slave[0] = shape_function_derivative_slave;

                    Geometry<Node<3>>::IntegrationPointsArrayType ips_slave(1);
                    ips_slave[0] = IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointType(local_coordinates[0], local_coordinates[1], integration_points[ip].weight);

                    IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointsContainerType IntegrationPointsArray_slave =
                    { { ips_slave } };
                    IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsValuesContainerType sfcontainer_slave =
                    { { shape_function_slave } };
                    IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsLocalGradientsContainerType sfgcontainer_slave =
                    { { shape_function_gradients_slave } };

                    Geometry<Node<3>>::Pointer intgeration_point_slave = Kratos::make_shared<IntegrationPointCurveOnSurface3d<Node<3>>>(
                        cps_slave,
                        IntegrationPointsArray_slave,
                        tangents_slave,
                        sfcontainer_slave,
                        sfgcontainer_slave);

                    Geometry<Node<3>>::Pointer coupling_geometry = Kratos::make_shared<CouplingMasterSlave<Node<3>>>(
                        intgeration_point_master, intgeration_point_slave, true);

                    new_geometries.push_back(coupling_geometry);
                }
            }
        }

        return new_geometries;
    }

    //static void GetIntegrationDomainCurve(
    //    ModelPart& rModelPart,
    //    std::shared_ptr<NodeCurveGeometry3D>& pNodeCurveGeometry3D,
    //    int& rShapeFunctionDerivativesOrder
    //    )
    //{
    //    auto spans = mNodeCurveGeometry3D->Spans();

    //    for (int i = 0; i < spans.size(); ++i)
    //    {
    //        ANurbs::Interval<double> domain(spans[i].T0(), spans[i].T1());

    //        int number_of_points = mNodeCurveGeometry3D->Degree() + 1;
    //        auto integration_points = ANurbs::IntegrationPoints<double>::Points1(number_of_points, domain);

    //        ANurbs::CurveShapeEvaluator<double> shape(mNodeCurveGeometry3D->Degree(), rShapeFunctionDerivativesOrder);

    //        for (int j = 0; j < integration_points.size(); ++j)
    //        {
    //            Vector local_coordinates(1);
    //            local_coordinates[0] = integration_points[j].t;

    //            shape.Compute(mNodeCurveGeometry3D->Knots(), integration_points[j].t);

    //            CreateElement(shape, rShapeFunctionDerivativesOrder);
    //        }
    //    }
    //}

    //static void GetIntegrationDomainPointOnCurve(
    //    ModelPart& rModelPart,
    //    const std::shared_ptr<NodeCurveGeometry3D>& pNodeCurveGeometry3D,
    //    const double& local_parameter
    //    )
    //{
    //    ANurbs::CurveShapeEvaluator<double> shape(
    //        mNodeCurveGeometry3D->Degree(),
    //        rShapeFunctionDerivativesOrder);

    //    shape.Compute(mNodeCurveGeometry3D->Knots(), local_parameter);

    //    Vector N_0 = ZeroVector(shape.NbNonzeroPoles());
    //    Matrix N_1 = ZeroMatrix(shape.NbNonzeroPoles(), 1);

    //    Element::GeometryType::PointsArrayType non_zero_control_points;
    //    for (int m = shape.FirstNonzeroPole(); m < shape.LastNonzeroPole() + 1; ++m)
    //    {
    //        non_zero_control_points.push_back(mNodeCurveGeometry3D->GetNode(m));

    //        N_0(m) = shape(0, m);
    //        N_1(m, 0) = shape(1, m);
    //    }
    //    if (rType == "element")
    //    {
    //        int id = 0;
    //        if (rModelPart.GetRootModelPart().Elements().size() > 0)
    //            id = rModelPart.GetRootModelPart().Elements().back().Id() + 1;

    //        auto element = rModelPart.CreateNewElement(
    //            rName,
    //            id,
    //            non_zero_control_points);

    //        element->SetValue(SHAPE_FUNCTION_VALUES, N_0);
    //        element->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
    //        element->SetValue(INTEGRATION_WEIGHT, integration_points[j].weight);

    //        element->SetValue(LOCAL_PARAMETER, local_parameter);
    //    }
    //    if (rType == "condition")
    //    {
    //        int id = 0;
    //        if (rModelPart.GetRootModelPart().Conditions().size() > 0)
    //            int id = rModelPart.GetRootModelPart().Conditions().back().Id() + 1;
    //        auto condition = rModelPart.CreateNewCondition(
    //            rName,
    //            id,
    //            non_zero_control_points);

    //        condition->SetValue(SHAPE_FUNCTION_VALUES, N_0);
    //        condition->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, N_1);
    //        condition->SetValue(INTEGRATION_WEIGHT, integration_points[j].weight);

    //        element->SetValue(LOCAL_PARAMETER, local_parameter);
    //    }
    //}

}

} // namespace Kratos

#endif // IGA_INTEGRATION_UTILITIES_H_INCLUDED
