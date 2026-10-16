//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti

#include "custom_utilities/brep_clipper_utilities.h"
#include "mapping_application.h"

namespace Kratos
{
namespace
{

template<class TLoopType>
void FillClipperPathFromLoop(
    const TLoopType& rLoop,
    Clipper2Lib::Path64& rPath,
    const double Factor,
    const double TessellationTolerance)
{
    Clipper2Lib::Point64 previous_point{
        std::numeric_limits<std::int64_t>::min(),
        std::numeric_limits<std::int64_t>::min()};

    for (IndexType i_curve = 0; i_curve < rLoop.size(); ++i_curve) {
        CurveTessellation<PointerVector<Node>> curve_tessellation;
        auto curve_geometry = *(rLoop[i_curve].get());
        curve_tessellation.Tessellate(
            curve_geometry,
            TessellationTolerance,
            1,
            true);

        std::vector<double> curve_spans;
        curve_geometry.SpansLocalSpace(curve_spans, 0);
        auto tessellation = curve_tessellation.GetTessellation();

        const bool is_degenerate_closed_tessellation =
            tessellation.size() == 2 &&
            norm_2(
                std::get<1>(tessellation.front()) -
                std::get<1>(tessellation.back())) < 1e-12;
        if (is_degenerate_closed_tessellation) {
            KRATOS_ERROR_IF(curve_spans.size() < 2)
                << "Cannot sample a closed trimming curve without a valid "
                << "parameter interval.\n";

            constexpr IndexType number_of_closed_curve_segments = 64;
            const double parameter_min = curve_spans.front();
            const double parameter_max = curve_spans.back();
            tessellation.clear();
            tessellation.reserve(number_of_closed_curve_segments + 1);

            for (IndexType i = 0;
                i <= number_of_closed_curve_segments;
                ++i) {
                const double position = static_cast<double>(i) /
                    static_cast<double>(number_of_closed_curve_segments);
                const double parameter =
                    (1.0 - position) * parameter_min +
                    position * parameter_max;

                array_1d<double, 3> surface_parameter = ZeroVector(3);
                surface_parameter[0] = parameter;
                curve_geometry.Calculate(
                    PARAMETER_2D_COORDINATES, surface_parameter);
                tessellation.emplace_back(parameter, surface_parameter);
            }
        }

        for (IndexType i_point = 0;
            i_point < tessellation.size();
            ++i_point) {
            const auto& r_uv = std::get<1>(tessellation[i_point]);
            const auto point = BrepTrimmingUtilities<false>::ToIntPoint(
                r_uv[0], r_uv[1], Factor);

            if (point.x != previous_point.x || point.y != previous_point.y) {
                rPath.push_back(point);
                previous_point = point;
            }
        }
    }
}

} // unnamed namespace

namespace BrepClipperUtilities
{

Clipper2Lib::Paths64 CreateAllLoops(
    const BrepSurface<PointerVector<Node>, false, PointerVector<Point>>& rBrepSurface,
    const double Factor,
    const double TessellationTolerance)
{
    const auto& r_outer_loops = rBrepSurface.GetOuterLoops();
    const auto& r_inner_loops = rBrepSurface.GetInnerLoops();
    KRATOS_ERROR_IF(r_outer_loops.empty())
        << "BrepSurface has no outer trimming loop.\n";

    Clipper2Lib::Paths64 all_loops(1 + r_inner_loops.size());
    FillClipperPathFromLoop(
        r_outer_loops[0], all_loops[0], Factor, TessellationTolerance);
    for (std::size_t i = 0; i < r_inner_loops.size(); ++i) {
        FillClipperPathFromLoop(
            r_inner_loops[i], all_loops[i + 1],
            Factor, TessellationTolerance);
    }
    return all_loops;
}

void SplitOuterAndInnerPaths(
    const Clipper2Lib::Paths64& rAllLoops,
    Clipper2Lib::Paths64& rOuterPaths,
    Clipper2Lib::Paths64& rInnerPaths)
{
    rOuterPaths.clear();
    rInnerPaths.clear();
    if (!rAllLoops.empty() && !rAllLoops[0].empty()) {
        rOuterPaths.push_back(rAllLoops[0]);
    }
    for (std::size_t i = 1; i < rAllLoops.size(); ++i) {
        if (!rAllLoops[i].empty()) {
            rInnerPaths.push_back(rAllLoops[i]);
        }
    }
}

Clipper2Lib::Paths64 ClipPathsWithTrimmedDomain(
    const Clipper2Lib::Paths64& rSubjectPaths,
    const Clipper2Lib::Paths64& rOuterPaths,
    const Clipper2Lib::Paths64& rInnerPaths)
{
    if (rSubjectPaths.empty() || rOuterPaths.empty()) {
        return {};
    }

    Clipper2Lib::Paths64 outer_solution;
    Clipper2Lib::Clipper64 outer_clipper;
    outer_clipper.AddSubject(rSubjectPaths);
    outer_clipper.AddClip(rOuterPaths);
    outer_clipper.Execute(
        Clipper2Lib::ClipType::Intersection,
        Clipper2Lib::FillRule::NonZero,
        outer_solution);

    if (outer_solution.empty() || rInnerPaths.empty()) {
        return outer_solution;
    }

    Clipper2Lib::Paths64 solution;
    Clipper2Lib::Clipper64 inner_clipper;
    inner_clipper.AddSubject(outer_solution);
    inner_clipper.AddClip(rInnerPaths);
    inner_clipper.Execute(
        Clipper2Lib::ClipType::Difference,
        Clipper2Lib::FillRule::NonZero,
        solution);
    return solution;
}

} // namespace BrepClipperUtilities
} // namespace Kratos
