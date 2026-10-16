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
#include <array>
#include <limits>
#include <string>
#include <vector>

// External includes

// Project includes
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "custom_utilities/volume_fraction_utilities.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"

namespace Kratos
{

namespace
{

using IndexType = VolumeFractionGsbmProcess::IndexType;

/// Closest point on the skins: the skin it lies on and the point
struct ClosestSkinPoint
{
    IndexType Skin = 0;
    array_1d<double, 3> Coordinates = ZeroVector(3);
};

/**
 * @brief Finds the closest point to rPoint over all segments of the skins, by
 * brute force.
 *
 * @param rSkins The skins, made of 2-noded line conditions.
 * @param rPoint The point to project.
 * @return The closest point with its skin. On ties the first one found is
 * kept, so the outer skin wins over the inner one.
 *
 * TODO: exact projection onto NURBS skins (e.g. as SBM's ProjectToNurbsSkin),
 * if the sampled polyline limits the accuracy of d.
 */
ClosestSkinPoint FindClosestSkinPoint(
    const std::array<const ModelPart*, 2>& rSkins,
    const array_1d<double, 3>& rPoint)
{
    ClosestSkinPoint closest;
    double min_distance = std::numeric_limits<double>::max();

    // Closest point of every segment of both skins, keeping the nearest one
    for (IndexType skin = 0; skin < rSkins.size(); ++skin) {
        for (const auto& r_segment : rSkins[skin]->Conditions()) {
            // Projection onto the segment's line, clamped to its end points
            const Point candidate = NearestPointUtilities::LineNearestPoint(
                rPoint, r_segment.GetGeometry()
            );
            const double distance = norm_2(candidate.Coordinates() - rPoint);

            if (distance < min_distance) {
                min_distance = distance;  // strictly nearer only
                closest = {skin, candidate.Coordinates()};
            }
        }
    }

    return closest;
}

} // namespace

void VolumeFractionGsbmProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const std::array<std::string, 4> required_parameters{
        "model_part_name", "geometry_name",
        "elements_model_part_name", "skeleton_model_part_name"
    };
    for (const std::string& r_key : required_parameters) {
        KRATOS_ERROR_IF(mParameters[r_key].GetString().empty())
            << "VolumeFractionGsbmProcess: \"" << r_key << "\" must be provided." << std::endl;
    }

    // Get the background patch
    const ModelPart& r_model_part = mpModel->GetModelPart(
        mParameters["model_part_name"].GetString()
    );
    const auto& r_background = r_model_part.GetGeometry(
        mParameters["geometry_name"].GetString()
    );
    KRATOS_ERROR_IF(r_background.LocalSpaceDimension() != 2)
        << "VolumeFractionGsbmProcess: only 2D background patches are supported." << std::endl;

    // Get the knot span boundaries, built such that uv = xy
    std::array<std::vector<double>, 2> spans;
    for (IndexType dir = 0; dir < 2; ++dir) {
        r_background.SpansLocalSpace(spans[dir], dir);
    }

    // Get the skins of the domain
    // TODO: NURBS skins, sampled here into line conditions (e.g. as
    // SnakeSbmProcess does); only polyline skins are supported for now
    Model empty_model;
    const ModelPart& r_empty_skin = empty_model.CreateModelPart("EmptySkin");
    const std::string outer_skin_name =
        mParameters["skin_model_part_outer_initial_name"].GetString();
    const std::string inner_skin_name =
        mParameters["skin_model_part_inner_initial_name"].GetString();
    const ModelPart& r_outer_skin = outer_skin_name.empty()
        ? r_empty_skin
        : mpModel->GetModelPart(outer_skin_name);
    const ModelPart& r_inner_skin = inner_skin_name.empty()
        ? r_empty_skin
        : mpModel->GetModelPart(inner_skin_name);

    // Get the volume fraction of every knot span
    const Matrix volume_fractions = VolumeFractionUtilities::ComputeKnotSpanVolumeFractions(
        r_outer_skin, r_inner_skin, spans[0], spans[1]
    );

    // Lambda for volume fraction assignment
    const auto assign_volume_fraction = [&volume_fractions](GeometryType& rGeometry) {
        KRATOS_ERROR_IF_NOT(rGeometry.Has(KNOT_SPAN_INDICES))
            << "VolumeFractionGsbmProcess: an entity has no KNOT_SPAN_INDICES. Elements and "
            << "facets must be the ones created by IgaModelerGsbm." << std::endl;
        const auto& r_indices = rGeometry.GetValue(KNOT_SPAN_INDICES);
        rGeometry.SetValue(VOLUME_FRACTION, volume_fractions(r_indices[0], r_indices[1]));
    };

    // Get the elements and assign volume fractions
    ModelPart& r_elements = mpModel->GetModelPart(
        mParameters["elements_model_part_name"].GetString()
    );
    for (auto& r_element : r_elements.Elements()) {
        assign_volume_fraction(r_element.GetGeometry());
    }

    // Get the facets and assign volume fractions
    ModelPart& r_facets = mpModel->GetModelPart(
        mParameters["skeleton_model_part_name"].GetString()
    );
    for (auto& r_facet : r_facets.Conditions()) {
        for (IndexType side = 0; side < 2; ++side) {
            // Part 0: - side (lower span); part 1: + side (upper span)
            assign_volume_fraction(r_facet.GetGeometry().GetGeometryPart(side));
        }
    }

    // Per skin (0: outer, 1: inner): the skin, its name and its identifier
    const std::array<const ModelPart*, 2> skins{&r_outer_skin, &r_inner_skin};
    const std::array<std::string, 2> skin_names{outer_skin_name, inner_skin_name};
    const std::array<std::string, 2> skin_identifiers{"outer", "inner"};

    // Per skin: where its projection nodes go, and the first free node Id
    std::array<ModelPart*, 2> projection_model_parts{nullptr, nullptr};
    IndexType node_id = 1;

    for (IndexType skin = 0; skin < 2; ++skin) {
        if (skin_names[skin].empty()) {
            continue;  // no skin of that kind, so nothing projects onto it
        }

        // Nodes at the closest skin points go to <skin>.projection_data, as SBM
        ModelPart& r_skin = mpModel->GetModelPart(skin_names[skin]);
        projection_model_parts[skin] = r_skin.HasSubModelPart("projection_data")
            ? &r_skin.GetSubModelPart("projection_data")
            : &r_skin.CreateSubModelPart("projection_data");

        // Node Ids are unique per root model part, so start after its largest
        const ModelPart& r_root_model_part = r_skin.GetRootModelPart();
        if (r_root_model_part.NumberOfNodes() > 0) {
            node_id = std::max(node_id, r_root_model_part.Nodes().back().Id() + 1);
        }
    }

    // Cut facets: points where alpha jumps, which get shifted boundary terms
    SizeType number_of_cut_facets = 0;

    for (auto& r_facet : r_facets.Conditions()) {
        const auto& r_minus = r_facet.GetGeometry().GetGeometryPart(0);
        const auto& r_plus = r_facet.GetGeometry().GetGeometryPart(1);

        if (r_minus.GetValue(VOLUME_FRACTION) == r_plus.GetValue(VOLUME_FRACTION)) {
            continue;  // same volume fraction on both sides: not cut, so no projection
        }
        ++number_of_cut_facets;

        // Closest skin point to the facet point, the same for both sides
        const ClosestSkinPoint closest = FindClosestSkinPoint(
            skins, r_minus.Center().Coordinates()
        );
        const IndexType skin = closest.Skin;

        // Node at the closest skin point, in that skin's projection data
        auto p_projection_node = projection_model_parts[skin]->CreateNewNode(
            node_id++, closest.Coordinates[0], closest.Coordinates[1], closest.Coordinates[2]
        );

        // The facet's link to its node, where the condition reads d and u_D
        r_facet.SetValue(NEIGHBOUR_NODES, GlobalPointersVector<Node>({p_projection_node}));
        // Which skin the node lies on
        r_facet.SetValue(IDENTIFIER, skin_identifiers[skin]);
    }

    // Number of cut knot spans with 0 < alpha < 1
    SizeType number_of_cut_spans = 0;
    for (IndexType i = 0; i < volume_fractions.size1(); ++i) {
        for (IndexType j = 0; j < volume_fractions.size2(); ++j) {
            const double alpha = volume_fractions(i, j);
            number_of_cut_spans += (alpha > 0.0 && alpha < 1.0);
        }
    }

    KRATOS_INFO_IF("VolumeFractionGsbmProcess", mParameters["echo_level"].GetInt() > 0)
        << "Assigned volume fractions to " << r_elements.NumberOfElements()
        << " elements and " << r_facets.NumberOfConditions() << " facets ("
        << number_of_cut_spans << " cut knot spans), and closest skin points to "
        << number_of_cut_facets << " cut facets." << std::endl;

    KRATOS_CATCH("")
}

const Parameters VolumeFractionGsbmProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "echo_level"                         : 0,
        "model_part_name"                    : "",
        "geometry_name"                      : "",
        "elements_model_part_name"           : "",
        "skeleton_model_part_name"           : "",
        "skin_model_part_outer_initial_name" : "",
        "skin_model_part_inner_initial_name" : ""
    })");
}

} // namespace Kratos
