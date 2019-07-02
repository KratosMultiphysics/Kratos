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

#if !defined(IGA_INTEGRATION_POINT_UTILITIES_H_INCLUDED)
#define IGA_INTEGRATION_POINT_UTILITIES_H_INCLUDED

#include "includes/model_part.h"
#include "geometries/integration_point_surface_3d.h"

namespace Kratos
{

namespace IgaIntegrationPointUtilities
{
    static Geometry<Node<3>>::Pointer GetIntegrationPointSurface(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
        const array_1d<double, 3>& rLocation,
        int ShapeFunctionDerivativesOrder)
    {
        std::vector<Geometry<Node<3>>::Pointer> new_geometries;

        auto spans_u = pSurface->SpansU();
        auto spans_v = pSurface->SpansV();

        ANurbs::SurfaceShapeEvaluator<double> shape(
            pSurface->DegreeU(),
            pSurface->DegreeV(),
            ShapeFunctionDerivativesOrder);

        shape.Compute(
            pSurface->KnotsU(),
            pSurface->KnotsV(),
            pSurface->Weights(),
            rLocation[0],
            rLocation[1]);

        Element::GeometryType::PointsArrayType control_points;
        Vector shape_function;
        Matrix shape_function_derivative;
        Matrix shape_function_second_derivative;

        Geometry<Node<3>>::PointsArrayType cps;

        int number_of_non_zero_cps = shape.NonzeroPoleIndices().size();

        Matrix N(1, number_of_non_zero_cps);
        Matrix DN_De(number_of_non_zero_cps, 2);

        array_1d<double, 3> location = ZeroVector(3);

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

        Geometry<Node<3>>::ShapeFunctionsGradientsType DN_De_type(1);
        DN_De_type[0] = DN_De;

        Geometry<Node<3>>::IntegrationPointsArrayType ips(1);
        ips[0] = IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointType(rLocation[0], rLocation[1], 1);

        IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointsContainerType ips_container =
        { { ips } };
        IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsValuesContainerType N_container =
        { { N } };
        IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsLocalGradientsContainerType DN_De_container =
        { { DN_De_type } };

        return std::make_shared<IntegrationPointSurface3d<Node<3>>>(
                cps,
                ips_container,
                N_container,
                DN_De_container);
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

#endif // IGA_INTEGRATION_POINT_UTILITIES_H_INCLUDED
