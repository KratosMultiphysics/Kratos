//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mario Caballero
//

// System includes
#include <algorithm>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <utility>

// External includes
#include "clipper/include/clipper2/clipper.h"

// Project includes
#include "custom_utilities/volume_fraction_utilities.h"

namespace Kratos
{

namespace
{

using IndexType = VolumeFractionUtilities::IndexType;

/// Maps physical coordinates to the integer lattice on which Clipper2 operates.
struct IntegerMapping
{
    double OriginX;
    double OriginY;
    double Scale;

    int64_t MapX(const double X) const
    {
        return static_cast<int64_t>(std::round((X - OriginX) * Scale));
    }

    int64_t MapY(const double Y) const
    {
        return static_cast<int64_t>(std::round((Y - OriginY) * Scale));
    }
};

/**
 * @brief Chains the 2-noded line conditions (e.g. LineCondition2D2N) of a skin model part into closed, positively oriented loops.
 *
 * @details Segments are chained by node connectivity, so neither the condition ordering nor the loop
 * orientation of the input matter.
 *
 * @param rLineConditions The 2-noded line conditions forming the skin.
 * @param rSkinModelPartName Name of the skin model part the conditions belong to, used in error messages.
 * @param rMapping Mapping from physical coordinates to Clipper2 integer units.
 */
Clipper2Lib::Paths64 ExtractClosedLoops(
    const ModelPart::ConditionsContainerType& rLineConditions,
    const std::string& rSkinModelPartName,
    const IntegerMapping& rMapping)
{
    // Map every segment start node Id to its node and to the Id of the segment end node
    std::unordered_map<IndexType, std::pair<const Node*, IndexType>> segments;

    for (const auto& r_condition : rLineConditions) {
        const auto& r_geometry = r_condition.GetGeometry();

        // Reject conditions that are not 2-noded lines
        KRATOS_ERROR_IF(r_geometry.size() != 2) << "Skin model part '" << rSkinModelPartName
            << "': condition " << r_condition.Id() << " has " << r_geometry.size()
            << " nodes. Only 2-noded line conditions are supported." << std::endl;

        // Fill the lookup table
        const bool is_new_start_node = segments.emplace(
            r_geometry[0].Id(), std::make_pair(&r_geometry[0], r_geometry[1].Id())).second;

        // Reject branching nodes: if this node already starts another segment
        KRATOS_ERROR_IF_NOT(is_new_start_node) << "Skin model part '" << rSkinModelPartName
            << "': node " << r_geometry[0].Id() << " starts more than one segment." << std::endl;
    }

    // Walk the lookup table loop by loop
    Clipper2Lib::Paths64 loops;
    // Remember which nodes are already part of a loop
    std::unordered_set<IndexType> visited_node_ids;

    for (const auto& r_condition : rLineConditions) {
        // Each condition is a candidate starting point
        const IndexType first_node_id = r_condition.GetGeometry()[0].Id();

        if (visited_node_ids.count(first_node_id) > 0)
            continue;  // skip early if a previous loop already went through it

        // Follow the segments node by node until the loop closes back at its first node
        Clipper2Lib::Path64 loop;
        IndexType current_node_id = first_node_id;

        do {
            const auto it_segment = segments.find(current_node_id);

            // Reject open chains: the boundary stops at a node where no segment starts
            KRATOS_ERROR_IF(it_segment == segments.end()) << "Skin model part '" << rSkinModelPartName
                << "': skin loop is not closed, no segment starts at node " << current_node_id << "." << std::endl;

            // Mark the node as used and append it to the loop in integer units
            visited_node_ids.insert(current_node_id);
            const Node& r_node = *(it_segment->second.first);
            loop.emplace_back(rMapping.MapX(r_node.X()), rMapping.MapY(r_node.Y()));

            // Move to the end node of the current segment
            current_node_id = it_segment->second.second;

            // Reject inner cycles: a node is reached again before the loop closes at its first node
            KRATOS_ERROR_IF(current_node_id != first_node_id && visited_node_ids.count(current_node_id) > 0)
                << "Skin model part '" << rSkinModelPartName << "': skin loop is not closed, node "
                << current_node_id << " is reached twice." << std::endl;
        } while (current_node_id != first_node_id);

        // Reject degenerate loops, which cannot enclose any area
        KRATOS_ERROR_IF(loop.size() < 3) << "Skin model part '" << rSkinModelPartName
            << "': skin loop starting at node " << first_node_id << " has less than 3 nodes." << std::endl;

        // Positive orientation for all loops, so that the NonZero fill rule does not depend on the input orientation
        if (Clipper2Lib::Area(loop) < 0.0)
            std::reverse(loop.begin(), loop.end());

        loops.push_back(std::move(loop));
    }

    return loops;
}

bool IsStrictlyIncreasing(const std::vector<double>& rValues)
{
    return std::adjacent_find(rValues.begin(), rValues.end(), std::greater_equal<double>()) == rValues.end();
}

} // namespace

Matrix VolumeFractionUtilities::ComputeKnotSpanVolumeFractions(
    const ModelPart& rOuterSkinModelPart,
    const ModelPart& rInnerSkinModelPart,
    const std::vector<double>& rSpansX,
    const std::vector<double>& rSpansY)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rSpansX.size() < 2 || rSpansY.size() < 2)
        << "At least two span boundaries are required in each direction." << std::endl;
    KRATOS_ERROR_IF_NOT(IsStrictlyIncreasing(rSpansX) && IsStrictlyIncreasing(rSpansY))
        << "Span boundaries must be strictly increasing." << std::endl;

    // Scale the integer factor to the largest extent of the span box
    const double box_extent = std::max(rSpansX.back() - rSpansX.front(), rSpansY.back() - rSpansY.front());
    const double factor = std::ldexp(1.0, 30) / box_extent;

    // Integer mapping for Clipper2
    const IntegerMapping mapping{rSpansX.front(), rSpansY.front(), factor};

    // Span boundaries in integer units
    std::vector<int64_t> spans_x(rSpansX.size());
    std::transform(rSpansX.begin(), rSpansX.end(), spans_x.begin(), [&mapping](const double X) { return mapping.MapX(X); });
    std::vector<int64_t> spans_y(rSpansY.size());
    std::transform(rSpansY.begin(), rSpansY.end(), spans_y.begin(), [&mapping](const double Y) { return mapping.MapY(Y); });

    // Turn the outer boundary into Clipper paths
    const Clipper2Lib::Paths64 outer_loops = rOuterSkinModelPart.NumberOfConditions() > 0
        // Outer skin given
        ? ExtractClosedLoops(rOuterSkinModelPart.Conditions(), rOuterSkinModelPart.FullName(), mapping)
        // No outer skin given, default to the rectangle
        : Clipper2Lib::Paths64{Clipper2Lib::Rect64(spans_x.front(), spans_y.front(), spans_x.back(), spans_y.back()).AsPath()};

    // Turn the inner boundary into Clipper paths
    const Clipper2Lib::Paths64 inner_loops = ExtractClosedLoops(rInnerSkinModelPart.Conditions(), rInnerSkinModelPart.FullName(), mapping);

    // Build the domain as the difference between the outer and inner loops
    const Clipper2Lib::Paths64 domain = Clipper2Lib::Difference(outer_loops, inner_loops, Clipper2Lib::FillRule::NonZero);
    // Build the bbox to skip spans that cannot intersect it
    const Clipper2Lib::Rect64 domain_bounds = Clipper2Lib::GetBounds(domain);

    // Number of non-zero knot spans (cells) in each direction
    const SizeType num_spans_x = rSpansX.size() - 1;
    const SizeType num_spans_y = rSpansY.size() - 1;

    // Volume fraction of span (i, j), zero-initialized for cells that never reach the domain
    Matrix volume_fractions = ZeroMatrix(num_spans_x, num_spans_y);

    for (IndexType i = 0; i < num_spans_x; ++i) {
        for (IndexType j = 0; j < num_spans_y; ++j) {
            // Build span (i,j)
            const Clipper2Lib::Rect64 span(spans_x[i], spans_y[j], spans_x[i + 1], spans_y[j + 1]);

            if (!span.Intersects(domain_bounds))
                continue;  // skip early if no intersection

            // Same rectangle as a polygon (its 4 corners), the input format of Clipper operations
            const Clipper2Lib::Path64 span_path = span.AsPath();
            // Part of the domain lying inside the span
            const Clipper2Lib::Paths64 cut = Clipper2Lib::Intersect({span_path}, domain, Clipper2Lib::FillRule::NonZero);

            // Volume fraction = cut area / span area
            // Holes have negative area, so they are subtracted
            volume_fractions(i, j) = std::clamp(Clipper2Lib::Area(cut) / Clipper2Lib::Area(span_path), 0.0, 1.0);
        }
    }

    return volume_fractions;

    KRATOS_CATCH("")
}

} // namespace Kratos
