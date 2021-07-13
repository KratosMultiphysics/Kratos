//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

#if !defined(KRATOS_BREP_TRIMMING_UTILITIES_H_INCLUDED)
#define KRATOS_BREP_TRIMMING_UTILITIES_H_INCLUDED

// System includes
#include "includes/define.h"

// External includes
#include "clipper/clipper.h"

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"

#include "utilities/tessellation_utilities/curve_tessellation.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    //using namespace ClipperLib;

    class BrepTrimmingUtilities
    {
    public:
        ///@name Type Definitions
        ///@{

        typedef std::size_t IndexType;
        typedef std::size_t SizeType;

        typedef IntegrationPoint<3> IntegrationPointType;
        typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

        typedef array_1d<double, 3> CoordinatesArrayType;

        //typedef CurveTessellation<PointerVector<Node<3>>> CurveTesselationType;
        typedef std::vector<std::pair<double, CoordinatesArrayType>> TessellationType;

        typedef signed long long cInt;

        template<class TBrepLoopType, class TPointType>
        static void CreateBrepSurfaceTrimmingIntegrationPoints(
            IntegrationPointsArrayType& rIntegrationPoints,
            const TBrepLoopType& rOuterLoops,
            const TBrepLoopType& rInnerLoops,
            const std::vector<double>& rSpansU,
            const std::vector<double>& rSpansV,
            IntegrationInfo& rIntegrationInfo,
            IndexType ClipperClipType = 0)
        {
            ClipperLib::Paths outer_loops(rOuterLoops.size()), inner_loops(rInnerLoops.size()), solution;
            const double factor = 1e-10;

            for (IndexType i = 0; i < rOuterLoops.size(); ++i) {
                for (IndexType j = 0; j < rOuterLoops[i].size(); ++j) {
                    CurveTessellation<PointerVector<TPointType>> curve_tesselation;
                    auto geometry_outer = *(rOuterLoops[i][j].get());
                    curve_tesselation.Tessellate(
                        geometry_outer, 0.01, 1, true);
                    auto tesselation = curve_tesselation.GetTessellation();
                    for (IndexType u = 0; u < tesselation.size(); ++u) {
                        outer_loops[i] << BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                    }
                }
            }

            for (IndexType i = 0; i < rInnerLoops.size(); ++i) {
                for (IndexType j = 0; j < rInnerLoops[i].size(); ++j) {
                    CurveTessellation<PointerVector<TPointType>> curve_tesselation;
                    auto geometry_inner = *(rInnerLoops[i][j].get());
                    curve_tesselation.Tessellate(
                        geometry_inner, 0.01, 1);
                    auto tesselation = curve_tesselation.GetTessellation();
                    for (IndexType u = 0; u < tesselation.size(); ++u) {
                        inner_loops[i] << BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                    }
                }
            }

            //perform intersection
            ClipperLib::Clipper c;
            ClipperLib::ClipType clipper_type = (ClipperClipType == 0)
                ? ClipperLib::ctIntersection
                : ClipperLib::ctDifference;
            ClipperLib::PolyType knots_polygon_type = (ClipperClipType == 0)
                ? ClipperLib::ptSubject
                : ClipperLib::ptClip;
            ClipperLib::PolyType trim_polygon_type = (ClipperClipType == 0)
                ? ClipperLib::ptClip
                : ClipperLib::ptSubject;

            for (IndexType i = 0; i < rSpansU.size() - 1; ++i) {
                for (IndexType j = 0; j < rSpansV.size() - 1; ++j) {

                    c.AddPaths(outer_loops, knots_polygon_type, true);
                    c.AddPaths(inner_loops, knots_polygon_type, true);

                    ClipperLib::Paths span(1);
                    span[0] <<
                        BrepTrimmingUtilities::ToIntPoint(rSpansU[i], rSpansV[j], factor) << BrepTrimmingUtilities::ToIntPoint(rSpansU[i + 1], rSpansV[j], factor) <<
                        BrepTrimmingUtilities::ToIntPoint(rSpansU[i + 1], rSpansV[j + 1], factor) << BrepTrimmingUtilities::ToIntPoint(rSpansU[i], rSpansV[j + 1], factor);

                    c.AddPaths(span, trim_polygon_type, true);
                    c.Execute(clipper_type, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

                    if (solution.size() == 0) {
                        continue;
                    }
                    else if (std::abs((std::abs(ClipperLib::Area(solution[0])) - std::abs(ClipperLib::Area(span[0]))) * (factor * factor)) < 1e-6) {

                        const IndexType number_of_integration_points = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) * rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1);

                        IndexType initial_integration_size = rIntegrationPoints.size();

                        if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                            rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                        }

                        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                        advance(integration_point_iterator, initial_integration_size);

                        IntegrationPointUtilities::IntegrationPoints2D(
                            integration_point_iterator,
                            rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0), rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1),
                            rSpansU[i], rSpansU[i + 1],
                            rSpansV[j], rSpansV[j + 1]);
                    }
                    else {
                        std::vector<Matrix> triangles;
                        BrepTrimmingUtilities::Triangulate_OPT(solution[0], triangles, factor);

                        const SizeType number_of_points = std::max(rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0), rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1));

                        const IndexType number_of_integration_points = triangles.size() * IntegrationPointUtilities::s_gauss_triangle[number_of_points].size();

                        IndexType initial_integration_size = rIntegrationPoints.size();

                        if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                            rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                        }

                        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                        advance(integration_point_iterator, initial_integration_size);

                        for (IndexType i = 0; i < triangles.size(); ++i)
                        {
                            IntegrationPointUtilities::IntegrationPointsTriangle2D(
                                integration_point_iterator,
                                number_of_points,
                                triangles[i](0, 0), triangles[i](1, 0), triangles[i](2, 0),
                                triangles[i](0, 1), triangles[i](1, 1), triangles[i](2, 1));
                        }
                    }
                    c.Clear();
                }
            }
        }

        template<class TBrepLoopType, class TPointType>
        static void CreateBrepSurfaceTrimmingIntegrationPointsReverseLoop(
            IntegrationPointsArrayType& rIntegrationPoints,
            const TBrepLoopType& rOuterLoops,
            const TBrepLoopType& rInnerLoops,
            const std::vector<double>& rSpansU,
            const std::vector<double>& rSpansV,
            IntegrationInfo& rIntegrationInfo)
        {
            ClipperLib::Paths outer_loops(rOuterLoops.size()), inner_loops(rInnerLoops.size()), solution;
            const double factor = 1e-10;

            for (IndexType i = 0; i < rOuterLoops.size(); ++i) {
                for (IndexType j = 0; j < rOuterLoops[i].size(); ++j) {
                    CurveTessellation<PointerVector<TPointType>> curve_tesselation;
                    auto geometry_outer = *(rOuterLoops[i][j].get());
                    curve_tesselation.Tessellate(
                        geometry_outer, 0.01, 1, true);
                    auto tesselation = curve_tesselation.GetTessellation();
                    for (
                        IndexType u = 0; u < tesselation.size(); ++u) {
                        outer_loops[i] << BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                    }
                }
            }

            for (IndexType i = 0; i < rInnerLoops.size(); ++i) {
                for (IndexType j = 0; j < rInnerLoops[i].size(); ++j) {
                    CurveTessellation<PointerVector<TPointType>> curve_tesselation;
                    auto geometry_inner = *(rInnerLoops[i][j].get());
                    curve_tesselation.Tessellate(
                        geometry_inner, 0.01, 1);
                    auto tesselation = curve_tesselation.GetTessellation();
                    for (IndexType u = 0; u < tesselation.size(); ++u) {
                        inner_loops[i] << BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                    }
                }
            }

            //perform intersection
            ClipperLib::Clipper c;

            //c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
            for (IndexType i = 0; i < rSpansU.size() - 1; ++i) {
                for (IndexType j = 0; j < rSpansV.size() - 1; ++j) {

                    c.AddPaths(outer_loops, ClipperLib::ptClip, true);
                    c.AddPaths(inner_loops, ClipperLib::ptClip, true);

                    ClipperLib::Paths span(1);
                    span[0] <<
                        BrepTrimmingUtilities::ToIntPoint(rSpansU[i], rSpansV[j], factor) << BrepTrimmingUtilities::ToIntPoint(rSpansU[i + 1], rSpansV[j], factor) <<
                        BrepTrimmingUtilities::ToIntPoint(rSpansU[i + 1], rSpansV[j + 1], factor) << BrepTrimmingUtilities::ToIntPoint(rSpansU[i], rSpansV[j + 1], factor);

                    c.AddPaths(span, ClipperLib::ptSubject, true);
                    c.Execute(ClipperLib::ctDifference, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

                    if (solution.size() == 0) {
                        continue;
                    }
                    else if (std::abs((std::abs(ClipperLib::Area(solution[0])) - std::abs(ClipperLib::Area(span[0]))) * (factor * factor)) < 1e-6) {
                        const IndexType number_of_integration_points = rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0) * rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1);

                        IndexType initial_integration_size = rIntegrationPoints.size();

                        if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                            rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                        }

                        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                        advance(integration_point_iterator, initial_integration_size);

                        IntegrationPointUtilities::IntegrationPoints2D(
                            integration_point_iterator,
                            rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0), rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1),
                            rSpansU[i], rSpansU[i + 1],
                            rSpansV[j], rSpansV[j + 1]);
                    }
                    else {
                        std::vector<Matrix> triangles;
                        BrepTrimmingUtilities::Triangulate_OPT(solution[0], triangles, factor);

                        const SizeType number_of_points = std::max(rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(0), rIntegrationInfo.GetNumberOfIntegrationPointsPerSpan(1));

                        const IndexType number_of_integration_points = triangles.size() * IntegrationPointUtilities::s_gauss_triangle[number_of_points].size();

                        IndexType initial_integration_size = rIntegrationPoints.size();

                        if (rIntegrationPoints.size() != initial_integration_size + number_of_integration_points) {
                            rIntegrationPoints.resize(initial_integration_size + number_of_integration_points);
                        }

                        typename IntegrationPointsArrayType::iterator integration_point_iterator = rIntegrationPoints.begin();
                        advance(integration_point_iterator, initial_integration_size);

                        for (IndexType i = 0; i < triangles.size(); ++i)
                        {
                            IntegrationPointUtilities::IntegrationPointsTriangle2D(
                                integration_point_iterator,
                                number_of_points,
                                triangles[i](0, 0), triangles[i](1, 0), triangles[i](2, 0),
                                triangles[i](0, 1), triangles[i](1, 1), triangles[i](2, 1));
                        }
                    }
                    c.Clear();
                }
            }
        }


        struct DPState {
            bool visible;
            double weight;
            long bestvertex;
        };
        struct Diagonal {
            long index1;
            long index2;
        };

        //Triangulation 
        static void Triangulate_OPT(const ClipperLib::Path& polygon, std::vector<Matrix>& triangles, const double factor)
        {
            array_1d<double, 2> p1, p2, p3, p4;
            int bestvertex;
            double weight = 0;
            double d1, d2 = 0.0;
            double minweight = std::numeric_limits<double>::max();
            Diagonal diagonal, newdiagonal;

            std::list<Diagonal> diagonals;

            IndexType n = polygon.size();
            std::vector< ClipperLib::IntPoint > const& points = polygon;
            matrix<DPState> dpstates(n, n);

            //init states and visibility
            for (IndexType i = 0; i < (n - 1); i++) {
                p1 = BrepTrimmingUtilities::IsConvex(points[i], factor);
                for (IndexType j = i + 1; j < n; j++) {
                    dpstates(j, i).visible = true;
                    dpstates(j, i).weight = 0;
                    dpstates(j, i).bestvertex = -1;
                    if (j != (i + 1)) {
                        p2 = BrepTrimmingUtilities::IsConvex(points[j], factor);

                        //visibility check
                        if (i == 0) p3 = BrepTrimmingUtilities::IsConvex(points[n - 1], factor);
                        else p3 = BrepTrimmingUtilities::IsConvex(points[i - 1], factor);
                        if (i == (n - 1)) p4 = BrepTrimmingUtilities::IsConvex(points[0], factor);
                        else p4 = BrepTrimmingUtilities::IsConvex(points[i + 1], factor);
                        if (!BrepTrimmingUtilities::InCone(p3, p1, p4, p2)) {
                            dpstates(j, i).visible = false;
                            continue;
                        }

                        if (j == 0) p3 = BrepTrimmingUtilities::IsConvex(points[n - 1], factor);
                        else p3 = BrepTrimmingUtilities::IsConvex(points[j - 1], factor);
                        if (j == (n - 1)) p4 = BrepTrimmingUtilities::IsConvex(points[0], factor);
                        else p4 = BrepTrimmingUtilities::IsConvex(points[j + 1], factor);
                        if (!BrepTrimmingUtilities::InCone(p3, p2, p4, p1)) {
                            dpstates(j, i).visible = false;
                            continue;
                        }

                        for (IndexType k = 0; k < n; k++) {
                            p3 = BrepTrimmingUtilities::IsConvex(points[k], factor);
                            if (k == (n - 1)) p4 = BrepTrimmingUtilities::IsConvex(points[0], factor);
                            else p4 = BrepTrimmingUtilities::IsConvex(points[k + 1], factor);
                            if (BrepTrimmingUtilities::Intersects(p1, p2, p3, p4)) {
                                dpstates(j, i).visible = false;
                                break;
                            }
                        }
                    }
                }
            }

            dpstates(n - 1, 0).visible = true;
            dpstates(n - 1, 0).weight = 0;
            dpstates(n - 1, 0).bestvertex = -1;

            for (IndexType gap = 2; gap < n; gap++) {
                for (IndexType i = 0; i < (n - gap); i++) {
                    IndexType j = i + gap;
                    if (!dpstates(j, i).visible) continue;
                    bestvertex = -1;
                    for (IndexType k = (i + 1); k < j; k++) {
                        if (!dpstates(k, i).visible) continue;
                        if (!dpstates(j, k).visible) continue;

                        if (k <= (i + 1)) d1 = 0;
                        else d1 = BrepTrimmingUtilities::Distance(BrepTrimmingUtilities::IsConvex(points[i], factor), BrepTrimmingUtilities::IsConvex(points[k], factor));
                        if (j <= (k + 1)) d2 = 0;
                        else d2 = BrepTrimmingUtilities::Distance(BrepTrimmingUtilities::IsConvex(points[k], factor), BrepTrimmingUtilities::IsConvex(points[j], factor));

                        weight = dpstates(k, i).weight + dpstates(j, k).weight + d1 + d2;

                        if ((bestvertex == -1) || (weight < minweight)) {
                            bestvertex = k;
                            minweight = weight;
                        }
                    }
                    if (bestvertex == -1) {
                        KRATOS_THROW_ERROR(std::runtime_error, "Triangulate: No points in polygon.", std::endl);
                    }

                    dpstates(j, i).bestvertex = bestvertex;
                    dpstates(j, i).weight = minweight;
                }
            }

            newdiagonal.index1 = 0;
            newdiagonal.index2 = n - 1;
            diagonals.push_back(newdiagonal);

            while (!diagonals.empty()) {
                diagonal = *(diagonals.begin());
                diagonals.pop_front();

                bestvertex = dpstates(diagonal.index2, diagonal.index1).bestvertex;
                if (bestvertex == -1) {
                    break;
                }
                Matrix triangle(3, 2);
                triangle(0, 0) = BrepTrimmingUtilities::IsConvex(points[diagonal.index1], factor)[0];
                triangle(0, 1) = BrepTrimmingUtilities::IsConvex(points[diagonal.index1], factor)[1];
                triangle(1, 0) = BrepTrimmingUtilities::IsConvex(points[bestvertex], factor)[0];
                triangle(1, 1) = BrepTrimmingUtilities::IsConvex(points[bestvertex], factor)[1];
                triangle(2, 0) = BrepTrimmingUtilities::IsConvex(points[diagonal.index2], factor)[0];
                triangle(2, 1) = BrepTrimmingUtilities::IsConvex(points[diagonal.index2], factor)[1];

                if (BrepTrimmingUtilities::GetAreaOfTriangle(triangle) > 1e-5)
                    triangles.push_back(triangle);
                else
                {
                    std::cout << "triangle with zero area" << BrepTrimmingUtilities::GetAreaOfTriangle(triangle) << std::endl;
                }
                if (bestvertex > (diagonal.index1 + 1)) {
                    newdiagonal.index1 = diagonal.index1;
                    newdiagonal.index2 = bestvertex;
                    diagonals.push_back(newdiagonal);
                }
                if (diagonal.index2 > (bestvertex + 1)) {
                    newdiagonal.index1 = bestvertex;
                    newdiagonal.index2 = diagonal.index2;
                    diagonals.push_back(newdiagonal);
                }
            }
        }

        static bool InCone(array_1d<double, 2>& p1, array_1d<double, 2>& p2,
            array_1d<double, 2>& p3, array_1d<double, 2>& p)
        {
            if (BrepTrimmingUtilities::IsConvex(p1, p2, p3)) {
                if (!BrepTrimmingUtilities::IsConvex(p1, p2, p)) return false;
                if (!BrepTrimmingUtilities::IsConvex(p2, p3, p)) return false;
                return true;
            }
            else {
                if (BrepTrimmingUtilities::IsConvex(p1, p2, p)) return true;
                if (BrepTrimmingUtilities::IsConvex(p2, p3, p)) return true;
                return false;
            }
        }

        static bool IsConvex(
            const array_1d<double, 2>& p1,
            const array_1d<double, 2>& p2,
            const array_1d<double, 2>& p3)
        {
            double tmp;
            tmp = (p3[1] - p1[1]) * (p2[0] - p1[0]) - (p3[0] - p1[0]) * (p2[1] - p1[1]);
            if (tmp > 0) return true;
            else return false;
        }

        static double Distance(array_1d<double, 2> point_1, array_1d<double, 2> point_2)
        {
            return sqrt(point_1[0] * point_2[0] + point_1[1] * point_2[1]);
        }

        static double GetAreaOfTriangle(const Matrix& triangle)
        {
            double area = (triangle(0, 0) * (triangle(1, 1) - triangle(2, 1))
                + triangle(1, 0) * (triangle(2, 1) - triangle(0, 1))
                + triangle(2, 0) * (triangle(0, 1) - triangle(1, 1))) / 2;

            return area;
        }

        //checks if two lines intersect
        static bool Intersects(array_1d<double, 2>& p11, array_1d<double, 2>& p12,
            array_1d<double, 2>& p21, array_1d<double, 2>& p22)
        {
            if ((p11[0] == p21[0]) && (p11[1] == p21[1])) return false;
            if ((p11[0] == p22[0]) && (p11[1] == p22[1])) return false;
            if ((p12[0] == p21[0]) && (p12[1] == p21[1])) return false;
            if ((p12[0] == p22[0]) && (p12[1] == p22[1])) return false;

            array_1d<double, 2> v1ort, v2ort, v;
            double dot11, dot12, dot21, dot22;

            v1ort[0] = p12[1] - p11[1];
            v1ort[1] = p11[0] - p12[0];

            v2ort[0] = p22[1] - p21[1];
            v2ort[1] = p21[0] - p22[0];

            v[0] = p21[0] - p11[0];
            v[1] = p21[1] - p11[1];
            dot21 = v[0] * v1ort[0] + v[1] * v1ort[1];
            v[0] = p22[0] - p11[0];
            v[1] = p22[1] - p11[1];
            dot22 = v[0] * v1ort[0] + v[1] * v1ort[1];

            v[0] = p11[0] - p21[0];
            v[1] = p11[1] - p21[1];
            dot11 = v[0] * v2ort[0] + v[1] * v2ort[1];
            v[0] = p12[0] - p21[0];
            v[1] = p12[1] - p21[1];
            dot12 = v[0] * v2ort[0] + v[1] * v2ort[1];

            if (dot11 * dot12 > 0) return false;
            if (dot21 * dot22 > 0) return false;

            return true;
        }

        static ClipperLib::IntPoint ToIntPoint(
            const double x,
            const double y,
            const double factor)
        {
            ClipperLib::IntPoint intPoint;

            intPoint.X = static_cast<cInt>(x / factor);
            intPoint.Y = static_cast<cInt>(y / factor);

            return intPoint;
        }

        static array_1d<double, 2> IsConvex(
            const ClipperLib::IntPoint& intPoint,
            const double factor)
        {
            array_1d<double, 2> point;

            point[0] = intPoint.X * factor;
            point[1] = intPoint.Y * factor;

            return point;
        }

        ///@}
    };
    ///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_UTILITIES_H_INCLUDED defined
