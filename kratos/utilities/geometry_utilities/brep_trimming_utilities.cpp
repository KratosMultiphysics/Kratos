//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// header includes
#include "brep_trimming_utilities.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    void BrepTrimmingUtilities::CreateBrepSurfaceTrimmingIntegrationPoints(
        IntegrationPointsArrayType& rIntegrationPoints,
        const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rOuterLoops,
        const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rInnerLoops,
        const std::vector<double>& rSpansU,
        const std::vector<double>& rSpansV,
        IntegrationInfo& rIntegrationInfo)
    {
        for (IndexType i_outer_loops = 0; i_outer_loops < rOuterLoops.size(); ++i_outer_loops) {

            Clipper2Lib::Paths64 all_loops(1 + rInnerLoops.size()), solution, solution_inner;
            const double factor = 1e-10;

            Clipper2Lib::Point64 int_point;
            int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
            int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());
            for (IndexType j = 0; j < rOuterLoops[i_outer_loops].size(); ++j) {
                CurveTessellation<PointerVector<Node>> curve_tesselation;
                auto geometry_outer = *(rOuterLoops[i_outer_loops][j].get());
                curve_tesselation.Tessellate(
                    geometry_outer, 0.001, 1, true);
                auto tesselation = curve_tesselation.GetTessellation();
                for (IndexType u = 0; u < tesselation.size(); ++u) {
                    auto new_int_point = BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                    if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y)) {
                        all_loops[i_outer_loops].push_back(new_int_point);
                        int_point.x = new_int_point.x;
                        int_point.y = new_int_point.y;
                    }
                }
            }

            for (IndexType i_inner_loops = 0; i_inner_loops < rInnerLoops.size(); ++i_inner_loops) {
                //ClipperLib::IntPoint int_point;
                int_point.x = static_cast<cInt>(std::numeric_limits<int>::min());
                int_point.y = static_cast<cInt>(std::numeric_limits<int>::min());
                for (IndexType j = 0; j < rInnerLoops[i_inner_loops].size(); ++j) {
                    CurveTessellation<PointerVector<Node>> curve_tesselation;
                    auto geometry_inner = *(rInnerLoops[i_inner_loops][j].get());
                    curve_tesselation.Tessellate(
                        geometry_inner, 0.001, 1, true);
                    auto tesselation = curve_tesselation.GetTessellation();
                    for (IndexType u = 0; u < tesselation.size(); ++u) {
                        auto new_int_point = BrepTrimmingUtilities::ToIntPoint(std::get<1>(tesselation[u])[0], std::get<1>(tesselation[u])[1], factor);
                        if (!(int_point.x == new_int_point.x && int_point.y == new_int_point.y)) {
                            all_loops[i_inner_loops + 1].push_back(new_int_point);
                            int_point.x = new_int_point.x;
                            int_point.y = new_int_point.x;
                        }
                    }
                }
            }

            for (IndexType i = 0; i < rSpansU.size() - 1; ++i) {
                for (IndexType j = 0; j < rSpansV.size() - 1; ++j) {
                    Clipper2Lib::Clipper64 c;
                    c.AddSubject(all_loops);

                    Clipper2Lib::Rect64 rectangle = Clipper2Lib::Rect64(
                        static_cast<cInt>(rSpansU[i] / factor), static_cast<cInt>(rSpansV[j] / factor),
                        static_cast<cInt>(rSpansU[i + 1] / factor), static_cast<cInt>(rSpansV[j + 1] / factor));

                    solution = Clipper2Lib::RectClip(rectangle, all_loops);

                    const double span_area = std::abs(Clipper2Lib::Area(rectangle.AsPath()));
                    double clip_area = 0.0;
                    if (solution.size() > 0)
                    {
                        clip_area = std::abs(Clipper2Lib::Area(solution[0]));
                        for (IndexType k = 1; k < solution.size(); ++k) {
                            clip_area -= std::abs(Clipper2Lib::Area(solution[k]));
                        }
                    }

                    Clipper2Lib::Clipper64 d;
                    d.AddSubject(solution);
                    d.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, solution_inner);

                    if (solution.size() == 0 || clip_area == 0.0) {
                        continue;
                    }
                    else if (std::abs(clip_area- span_area) < 1000) {
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
                        for(IndexType i = 0; i < solution_inner.size(); ++i)
                        {
                            BrepTrimmingUtilities::Triangulate_OPT(solution_inner[i], triangles, factor);    
                        }
                        

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
    };
    ///@} // Kratos Classes

    //template void BrepTrimmingUtilities::CreateBrepSurfaceTrimmingIntegrationPoints<
    //    DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>, Node>(
    //    BrepTrimmingUtilities::IntegrationPointsArrayType& rIntegrationPoints,
    //    const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rOuterLoops,
    //    const DenseVector<DenseVector<typename BrepCurveOnSurface<PointerVector<Node>, PointerVector<Point>>::Pointer>>& rInnerLoops,
    //    const std::vector<double>& rSpansU,
    //    const std::vector<double>& rSpansV,
    //    IntegrationInfo& rIntegrationInfo);

} // namespace Kratos.
