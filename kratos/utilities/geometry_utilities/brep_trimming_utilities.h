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

// Std includes
#include <list>

// System includes
#include "includes/define.h"

// External includes
#include "clipper/include/clipper2/clipper.h"

// Project includes
#include "geometries/geometry.h"
#include "geometries/brep_curve_on_surface.h"

#include "utilities/tessellation_utilities/curve_tessellation.h"
#include "includes/node.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    //using namespace ClipperLib;

    class KRATOS_API(KRATOS_CORE) BrepTrimmingUtilities
    {
    public:
        ///@name Type Definitions
        ///@{

        typedef std::size_t IndexType;
        typedef std::size_t SizeType;

        typedef IntegrationPoint<3> IntegrationPointType;
        typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;

        typedef array_1d<double, 3> CoordinatesArrayType;

        //typedef CurveTessellation<PointerVector<Node>> CurveTesselationType;
        typedef std::vector<std::pair<double, CoordinatesArrayType>> TessellationType;

        typedef signed long long cInt;

        using BrepCurveOnSurfacePointerType = typename BrepCurveOnSurface<PointerVector<Node>, false, PointerVector<Point>>::Pointer;

        //template<class TBrepLoopType, class TPointType>
        static void CreateBrepSurfaceTrimmingIntegrationPoints(
            IntegrationPointsArrayType& rIntegrationPoints,
            const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rOuterLoops,
            const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rInnerLoops,
            const std::vector<double>& rSpansU,
            const std::vector<double>& rSpansV,
            IntegrationInfo& rIntegrationInfo);

        struct DPState {
            bool visible;
            double weight;
            long best_vertex;
        };

        struct Diagonal {
            long index1;
            long index2;
        };

        //Triangulation 
        static void Triangulate_OPT(const Clipper2Lib::Path64& polygon, std::vector<Matrix>& triangles, const double factor)
        {
            array_1d<double, 2> p1, p2, p3, p4;
            int best_vertex;
            double weight = 0;
            double d1 = 0.0, d2 = 0.0;
            double min_weight = std::numeric_limits<double>::max();
            Diagonal diagonal, new_diagonal;

            std::queue<Diagonal> diagonals;

            IndexType n = polygon.size();
            std::vector< Clipper2Lib::Point64 > const& points = polygon;

            //if first and last point are coincide, neglect the last point
            double p0_x, p0_y, pn_x, pn_y, dpx, dpy;
            p0_x = BrepTrimmingUtilities::IntPointToDoublePoint(points[0], factor)[0];
            p0_y = BrepTrimmingUtilities::IntPointToDoublePoint(points[0], factor)[1];
            pn_x = BrepTrimmingUtilities::IntPointToDoublePoint(points[n - 1], factor)[0];
            pn_y = BrepTrimmingUtilities::IntPointToDoublePoint(points[n - 1], factor)[1];
            dpx = pn_x - p0_x;
            dpy = pn_y - p0_y;

            if(sqrt((dpx*dpx+dpy*dpy)) < 1e-9){
                n = n - 1;
            }

            if (n == 3) //special case with only one triangle
            {
                Matrix triangle(3, 2);
                triangle(0, 0) = BrepTrimmingUtilities::IntPointToDoublePoint(polygon[0], factor)[0];
                triangle(0, 1) = BrepTrimmingUtilities::IntPointToDoublePoint(polygon[0], factor)[1];
                triangle(1, 0) = BrepTrimmingUtilities::IntPointToDoublePoint(polygon[1], factor)[0];
                triangle(1, 1) = BrepTrimmingUtilities::IntPointToDoublePoint(polygon[1], factor)[1];
                triangle(2, 0) = BrepTrimmingUtilities::IntPointToDoublePoint(polygon[2], factor)[0];
                triangle(2, 1) = BrepTrimmingUtilities::IntPointToDoublePoint(polygon[2], factor)[1];
                triangles.push_back(triangle);
                return;
            }

            matrix<DPState> dpstates(n, n);

            //init states and visibility
            for (IndexType i = 0; i < (n - 1); i++) {
                p1 = BrepTrimmingUtilities::IntPointToDoublePoint(points[i], factor);
                for (IndexType j = i + 1; j < n; j++) {
                    dpstates(j, i).visible = true;
                    dpstates(j, i).weight = 0;
                    dpstates(j, i).best_vertex = -1;
                    if (j != (i + 1)) {
                        p2 = BrepTrimmingUtilities::IntPointToDoublePoint(points[j], factor);

                        //visibility check
                        if (i == 0) p3 = BrepTrimmingUtilities::IntPointToDoublePoint(points[n - 1], factor);
                        else p3 = BrepTrimmingUtilities::IntPointToDoublePoint(points[i - 1], factor);
                        if (i == (n - 1)) p4 = BrepTrimmingUtilities::IntPointToDoublePoint(points[0], factor);
                        else p4 = BrepTrimmingUtilities::IntPointToDoublePoint(points[i + 1], factor);
                        if (!BrepTrimmingUtilities::InCone(p3, p1, p4, p2)) {
                            dpstates(j, i).visible = false;
                            continue;
                        }

                        if (j == 0) p3 = BrepTrimmingUtilities::IntPointToDoublePoint(points[n - 1], factor);
                        else p3 = BrepTrimmingUtilities::IntPointToDoublePoint(points[j - 1], factor);
                        if (j == (n - 1)) p4 = BrepTrimmingUtilities::IntPointToDoublePoint(points[0], factor);
                        else p4 = BrepTrimmingUtilities::IntPointToDoublePoint(points[j + 1], factor);
                        if (!BrepTrimmingUtilities::InCone(p3, p2, p4, p1)) {
                            dpstates(j, i).visible = false;
                            continue;
                        }

                        for (IndexType k = 0; k < n; k++) {
                            p3 = BrepTrimmingUtilities::IntPointToDoublePoint(points[k], factor);
                            if (k == (n - 1)) p4 = BrepTrimmingUtilities::IntPointToDoublePoint(points[0], factor);
                            else p4 = BrepTrimmingUtilities::IntPointToDoublePoint(points[k + 1], factor);
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
            dpstates(n - 1, 0).best_vertex = -1;

            for (IndexType gap = 2; gap < n; gap++) {
                for (IndexType i = 0; i < (n - gap); i++) {
                    IndexType j = i + gap;
                    if (!dpstates(j, i).visible) continue;
                    best_vertex = -1;
                    for (IndexType k = (i + 1); k < j; k++) {
                        if (!dpstates(k, i).visible) continue;
                        if (!dpstates(j, k).visible) continue;

                        if (k <= (i + 1)) d1 = 0;
                        else d1 = BrepTrimmingUtilities::Distance(BrepTrimmingUtilities::IntPointToDoublePoint(points[i], factor), BrepTrimmingUtilities::IntPointToDoublePoint(points[k], factor));
                        if (j <= (k + 1)) d2 = 0;
                        else d2 = BrepTrimmingUtilities::Distance(BrepTrimmingUtilities::IntPointToDoublePoint(points[k], factor), BrepTrimmingUtilities::IntPointToDoublePoint(points[j], factor));

                        weight = dpstates(k, i).weight + dpstates(j, k).weight + d1 + d2;

                        if ((best_vertex == -1) || (weight < min_weight)) {
                            best_vertex = k;
                            min_weight = weight;
                        }
                    }
                    if (best_vertex == -1) {
                        return;
                    }

                    dpstates(j, i).best_vertex = best_vertex;
                    dpstates(j, i).weight = min_weight;
                }
            }

            new_diagonal.index1 = 0;
            new_diagonal.index2 = n - 1;
            diagonals.push(new_diagonal);

            while (!diagonals.empty()) {
                diagonal = diagonals.front();
                diagonals.pop();

                best_vertex = dpstates(diagonal.index2, diagonal.index1).best_vertex;
                if (best_vertex == -1) {
                    break;
                }
                Matrix triangle(3, 2);
                triangle(0, 0) = BrepTrimmingUtilities::IntPointToDoublePoint(points[diagonal.index1], factor)[0];
                triangle(0, 1) = BrepTrimmingUtilities::IntPointToDoublePoint(points[diagonal.index1], factor)[1];
                triangle(1, 0) = BrepTrimmingUtilities::IntPointToDoublePoint(points[best_vertex], factor)[0];
                triangle(1, 1) = BrepTrimmingUtilities::IntPointToDoublePoint(points[best_vertex], factor)[1];
                triangle(2, 0) = BrepTrimmingUtilities::IntPointToDoublePoint(points[diagonal.index2], factor)[0];
                triangle(2, 1) = BrepTrimmingUtilities::IntPointToDoublePoint(points[diagonal.index2], factor)[1];

                if (BrepTrimmingUtilities::GetAreaOfTriangle(triangle) > 1e-5)
                    triangles.push_back(triangle);
                else
                {
                    std::cout << "triangle with zero area" << BrepTrimmingUtilities::GetAreaOfTriangle(triangle) << std::endl;
                }
                if (best_vertex > (diagonal.index1 + 1)) {
                    new_diagonal.index1 = diagonal.index1;
                    new_diagonal.index2 = best_vertex;
                    diagonals.push(new_diagonal);
                }
                if (diagonal.index2 > (best_vertex + 1)) {
                    new_diagonal.index1 = best_vertex;
                    new_diagonal.index2 = diagonal.index2;
                    diagonals.push(new_diagonal);
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

        static Clipper2Lib::Point64 ToIntPoint(
            const double x,
            const double y,
            const double factor)
        {
            Clipper2Lib::Point64 int_point;

            int_point.x = static_cast<cInt>(x / factor);
            int_point.y = static_cast<cInt>(y / factor);

            return int_point;
        }

        static array_1d<double, 2> IntPointToDoublePoint(
            const Clipper2Lib::Point64& int_point,
            const double factor)
        {
            array_1d<double, 2> point;

            point[0] = int_point.x * factor;
            point[1] = int_point.y * factor;

            return point;
        }

        ///@}
    };
    ///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_BREP_TRIMMING_UTILITIES_H_INCLUDED defined
