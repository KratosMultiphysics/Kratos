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
    static int GetDerivativesSize(const int ShapeFunctionDerivativesOrder)
    {
        int derivatives_size = 0;
        for (int i = 0; i < ShapeFunctionDerivativesOrder; ++i)
        {
            derivatives_size += i + 1;
        }
        return derivatives_size;
    }

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

        Geometry<Node<3>>::PointsArrayType cps;

        int number_of_non_zero_cps = shape.NonzeroPoleIndices().size();

        int derivatives_size = GetDerivativesSize(ShapeFunctionDerivativesOrder);

        Matrix N(derivatives_size, number_of_non_zero_cps);

        array_1d<double, 3> location = ZeroVector(3);

        for (int n = 0; n < number_of_non_zero_cps; ++n)
        {
            int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
            int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

            cps.push_back(pSurface->GetNode(
                shape.NonzeroPoleIndices()[n].first,
                shape.NonzeroPoleIndices()[n].second));

            for (int i = 0; i < derivatives_size; ++i)
            {
                N(i, n) = shape(i, indexU, indexV);
            }
        }

        auto ip = IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointType(rLocation[0], rLocation[1], 1);
        Geometry<Node<3>>::IntegrationPointsArrayType ips(1);
        ips[0] = ip;

        IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointsContainerType ips_container =
        { { ips } };

        auto ar = DenseVector<Matrix>(1);
        ar[0] = N;

        //IntegrationPointCurveOnSurface3d<Node<3>>::ShapeFunctionsContainerType N_container =
        //{{ar}};
        return std::make_shared<IntegrationPointSurface3d<Node<3>>>(
            cps);
    }

    static std::vector<Geometry<Node<3>>::Pointer> GetIntegrationDomainGeometrySurface(
        const std::shared_ptr<NodeSurfaceGeometry3D>& pSurface,
        const TrimmedSurfaceClipping& rClipper,
        int ShapeFunctionDerivativesOrder)
    {
        std::vector<Geometry<Node<3>>::Pointer> new_geometries;

        int degree_u = pSurface->DegreeU();
        int degree_v = pSurface->DegreeV();

        int degree = std::max(degree_u, degree_v) + 1;

        ANurbs::SurfaceShapeEvaluator<double> shape(
            pSurface->DegreeU(),
            pSurface->DegreeV(),
            ShapeFunctionDerivativesOrder);

        int derivatives_size = GetDerivativesSize(ShapeFunctionDerivativesOrder);

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
                        shape.Compute(
                            pSurface->KnotsU(),
                            pSurface->KnotsV(),
                            pSurface->Weights(),
                            integration_points[i].u,
                            integration_points[i].v);

                        int number_of_cps = shape.NonzeroPoleIndices().size();

                        Element::GeometryType::PointsArrayType cps;
                        cps.reserve(number_of_cps);
                        Matrix N(derivatives_size, number_of_cps);

                        for (int n = 0; n < number_of_cps; ++n)
                        {
                            int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                            int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                            cps.push_back(pSurface->GetNode(
                                shape.NonzeroPoleIndices()[n].first,
                                shape.NonzeroPoleIndices()[n].second));

                            for (int j = 0; j < derivatives_size; ++j)
                            {
                                N(j, n) = shape(j, indexU, indexV);
                            }
                        }

                        auto ip = IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointType(
                            integration_points[i].u,
                            integration_points[i].v,
                            integration_points[i].weight);

                        //new_geometries.push_back(std::make_shared<IntegrationPointSurface3d<Node<3>>>(
                        //    cps, ip, N));
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
                            shape.Compute(
                                pSurface->KnotsU(),
                                pSurface->KnotsV(),
                                pSurface->Weights(),
                                integration_point_polygon.IntegrationPoint(i).u,
                                integration_point_polygon.IntegrationPoint(i).v);

                            int number_of_cps = shape.NonzeroPoleIndices().size();

                            Element::GeometryType::PointsArrayType cps;
                            cps.reserve(number_of_cps);
                            Matrix N(derivatives_size, number_of_cps);

                            for (int n = 0; n < number_of_cps; ++n)
                            {
                                int indexU = shape.NonzeroPoleIndices()[n].first - shape.FirstNonzeroPoleU();
                                int indexV = shape.NonzeroPoleIndices()[n].second - shape.FirstNonzeroPoleV();

                                cps.push_back(pSurface->GetNode(
                                    shape.NonzeroPoleIndices()[n].first,
                                    shape.NonzeroPoleIndices()[n].second));

                                for (int j = 0; j < derivatives_size; ++j)
                                {
                                    N(j, n) = shape(j, indexU, indexV);
                                }
                            }

                            auto ip = IntegrationPointCurveOnSurface3d<Node<3>>::IntegrationPointType(
                                integration_point_polygon.IntegrationPoint(i).u,
                                integration_point_polygon.IntegrationPoint(i).v,
                                integration_point_polygon.IntegrationPoint(i).weight);

                            //new_geometries.push_back(std::make_shared<IntegrationPointSurface3d<Node<3>>>(
                            //    cps, ip, N));
                        }
                    }
                }
            }
        }
        return new_geometries;
    }
}

} // namespace Kratos

#endif // IGA_INTEGRATION_POINT_UTILITIES_H_INCLUDED
