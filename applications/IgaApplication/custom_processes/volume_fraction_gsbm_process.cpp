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
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// External includes

// Project includes
#include "custom_modelers/iga_modeler_sbm.h"
#include "custom_processes/volume_fraction_gsbm_process.h"
#include "custom_utilities/volume_fraction_utilities.h"
#include "geometries/nurbs_curve_geometry.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"
#include "utilities/tessellation_utilities/curve_tessellation.h"

namespace Kratos
{

namespace
{

using SizeType = VolumeFractionGsbmProcess::SizeType;
using IndexType = VolumeFractionGsbmProcess::IndexType;
using GeometryType = VolumeFractionGsbmProcess::GeometryType;
using NurbsProjectionResult = IgaModelerSbm::NurbsProjectionResult;
using NurbsCurveType = NurbsCurveGeometry<2, PointerVector<Node>>;

/**
 * @brief Tessellates the NURBS curves of a skin into 2-noded line conditions,
 * up to a chordal error. Curves sharing an endpoint share its node, so the
 * curves of each loop chain head to tail. The curves are stored on rPolyline
 * (NEIGHBOUR_GEOMETRIES) for the projection, as SnakeSbmProcess does.
 *
 * @param rSkin The skin, with one NurbsCurveGeometry per curve.
 * @param rPolyline Where the nodes and line conditions are created.
 * @param Tolerance Chordal error of the tessellation.
 */
void TessellateNurbsSkin(
    const ModelPart& rSkin,
    ModelPart& rPolyline,
    const double Tolerance)
{
    // Ids after the largest ones of the root model part
    const ModelPart& r_root_model_part = rPolyline.GetRootModelPart();
    IndexType node_id = r_root_model_part.NumberOfNodes() > 0
        ? r_root_model_part.Nodes().back().Id() + 1
        : 1;
    IndexType condition_id = r_root_model_part.NumberOfConditions() > 0
        ? r_root_model_part.Conditions().back().Id() + 1
        : 1;
    // Recursively: ImportNurbsSbmModeler already creates them in the skin
    const Properties::Pointer p_properties = rPolyline.RecursivelyHasProperties(0)
        ? rPolyline.pGetProperties(0)
        : rPolyline.CreateNewProperties(0);

    // Curve endpoints, one node per point where curves meet
    std::vector<Node::Pointer> endpoints;
    const auto get_endpoint = [&](const array_1d<double, 3>& rPoint) {
        for (const auto& p_endpoint : endpoints) {
            if (norm_2(p_endpoint->Coordinates() - rPoint) < 1.0e-7) {
                return p_endpoint;
            }
        }
        endpoints.push_back(
            rPolyline.CreateNewNode(node_id++, rPoint[0], rPoint[1], rPoint[2])
        );
        return endpoints.back();
    };

    std::vector<GeometryType::Pointer> nurbs_curves;
    for (const auto& r_geometry : rSkin.Geometries()) {
        const auto p_curve = dynamic_cast<const NurbsCurveType*>(&r_geometry);
        KRATOS_ERROR_IF(p_curve == nullptr) << "VolumeFractionGsbmProcess: geometry #"
            << r_geometry.Id() << " of skin " << rSkin.FullName()
            << " is not a 2D NURBS curve." << std::endl;
        KRATOS_ERROR_IF_NOT(r_geometry.Has(IDENTIFIER)) << "VolumeFractionGsbmProcess: "
            << "curve #" << r_geometry.Id() << " of skin " << rSkin.FullName()
            << " has no IDENTIFIER (its layer), which the projection needs." << std::endl;
        nurbs_curves.push_back(rSkin.pGetGeometry(r_geometry.Id()));

        // Points of the curve up to the chordal error, end points included
        const auto points = CurveTessellation<PointerVector<Node>>::ComputeTessellation(
            *p_curve, p_curve->PolynomialDegree(0), p_curve->DomainInterval(),
            p_curve->KnotSpanIntervals(), Tolerance
        );

        // Shared nodes at the endpoints, new ones inside
        std::vector<Node::Pointer> nodes{get_endpoint(points.front().second)};
        for (IndexType i = 1; i + 1 < points.size(); ++i) {
            const auto& r_point = points[i].second;
            nodes.push_back(
                rPolyline.CreateNewNode(node_id++, r_point[0], r_point[1], r_point[2])
            );
        }
        nodes.push_back(get_endpoint(points.back().second));

        // One segment per pair of consecutive nodes
        for (IndexType i = 0; i + 1 < nodes.size(); ++i) {
            const std::vector<IndexType> node_ids{nodes[i]->Id(), nodes[i + 1]->Id()};
            rPolyline.CreateNewCondition(
                "LineCondition2D2N", condition_id++, node_ids, p_properties
            );
        }
    }
    rPolyline.SetValue(NEIGHBOUR_GEOMETRIES, nurbs_curves);
}

/**
 * @brief Finds the closest point to rPoint on a skin: exactly on NURBS skins,
 * with SBM's ProjectToNurbsSkin, and by brute force over the segments of
 * polylines.
 *
 * @param rSkin The skin's polyline, of 2-noded line conditions.
 * @param rPoint The point to project.
 * @return The closest point, as SBM's projection result. For polylines only
 * pGeometry (the segment), SkinCoordinates and SquaredDistance are set; an
 * empty skin keeps the largest SquaredDistance.
 */
NurbsProjectionResult FindClosestSkinPoint(
    const ModelPart& rSkin,
    const array_1d<double, 3>& rPoint)
{
    // NURBS skins: Newton to round-off, over all curves (no search radius)
    if (rSkin.Has(NEIGHBOUR_GEOMETRIES)) {
        return IgaModelerSbm::ProjectToNurbsSkin(
            rPoint, rSkin.GetValue(NEIGHBOUR_GEOMETRIES),
            std::numeric_limits<double>::infinity(), 1.0e-12, false
        );
    }

    // Polylines: closest point of every segment, keeping the nearest one
    NurbsProjectionResult closest;
    for (const auto& r_segment : rSkin.Conditions()) {
        // Projection onto the segment's line, clamped to its end points
        const Point candidate = NearestPointUtilities::LineNearestPoint(
            rPoint, r_segment.GetGeometry()
        );
        const array_1d<double, 3> distance = candidate.Coordinates() - rPoint;
        const double squared_distance = inner_prod(distance, distance);

        if (squared_distance < closest.SquaredDistance) {  // strictly nearer only
            closest.pGeometry = r_segment.pGetGeometry();
            closest.SkinCoordinates = candidate.Coordinates();
            closest.SquaredDistance = squared_distance;
        }
    }

    return closest;
}

/**
 * @brief Stores the boundary curve's geometry at the closest point on rNode:
 * the outward NORMAL \nu, the LOCAL_TANGENT \tau with the domain on its left,
 * the CURVATURE \kappa and CURVATURE_DERIVATIVES [\kappa', \kappa''] in arc
 * length, with \tau' = -\kappa \nu.
 *
 * @param rNode The node at the closest point.
 * @param rProjection The closest point, with its curve or segment.
 * @param Sign +1 on the outer skin, -1 on the inner one, as SBM's normal_sign.
 */
void AssignBoundaryCurveGeometry(
    Node& rNode,
    const NurbsProjectionResult& rProjection,
    const double Sign)
{
    array_1d<double, 3> tangent;
    double curvature = 0.0;
    Vector curvature_derivatives = ZeroVector(2);

    const auto p_curve = dynamic_cast<const NurbsCurveType*>(rProjection.pGeometry.get());
    if (p_curve == nullptr) {
        const auto& r_line = *rProjection.pGeometry;
        tangent = Sign * (r_line[1].Coordinates() - r_line[0].Coordinates()) / r_line.Length();
    } else {
        // c, c', c'', c''' and c'''' of the curve at the closest point
        std::vector<array_1d<double, 3>> c;
        p_curve->GlobalSpaceDerivatives(c, rProjection.LocalCoordinates, 4);

        // A = c' x c'' and S = |c'|^2, with their u-derivatives
        const double a = c[1][0] * c[2][1] - c[1][1] * c[2][0];
        const double a_u = c[1][0] * c[3][1] - c[1][1] * c[3][0];
        const double a_uu = c[1][0] * c[4][1] - c[1][1] * c[4][0]
            + c[2][0] * c[3][1] - c[2][1] * c[3][0];
        const double s = inner_prod(c[1], c[1]);
        const double s_u = 2.0 * inner_prod(c[1], c[2]);
        const double s_uu = 2.0 * (inner_prod(c[2], c[2]) + inner_prod(c[1], c[3]));

        // \kappa_c = A / S^(3/2) along the curve, and its u-derivatives
        const double s_3_2 = s * std::sqrt(s);
        const double kappa = a / s_3_2;
        const double kappa_u = (a_u - 1.5 * a * s_u / s) / s_3_2;
        const double kappa_uu = (a_uu - 3.0 * a_u * s_u / s
            + 3.75 * a * s_u * s_u / (s * s) - 1.5 * a * s_uu / s) / s_3_2;

        // d/ds = S^(-1/2) d/du, and \kappa^(j) = Sign^(j+1) \kappa_c^(j)
        tangent = Sign * c[1] / std::sqrt(s);
        curvature = Sign * kappa;
        curvature_derivatives[0] = kappa_u / std::sqrt(s);
        curvature_derivatives[1] = Sign * (kappa_uu - 0.5 * kappa_u * s_u / s) / s;
    }

    // \nu = (\tau_y, -\tau_x): the right of \tau, out of the domain
    array_1d<double, 3> normal = ZeroVector(3);
    normal[0] = tangent[1];
    normal[1] = -tangent[0];

    rNode.SetValue(NORMAL, normal);
    rNode.SetValue(LOCAL_TANGENT, tangent);
    rNode.SetValue(CURVATURE, curvature);
    rNode.SetValue(CURVATURE_DERIVATIVES, curvature_derivatives);
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
    Model empty_model;
    const ModelPart& r_empty_skin = empty_model.CreateModelPart("EmptySkin");
    const std::string outer_skin_name =
        mParameters["skin_model_part_outer_initial_name"].GetString();
    const std::string inner_skin_name =
        mParameters["skin_model_part_inner_initial_name"].GetString();
    const double tolerance = mParameters["tessellation_tolerance"].GetDouble();

    // Polyline of a skin: itself, or <skin>.tessellation for NURBS curves
    const auto get_polyline = [&](const std::string& rName) -> const ModelPart& {
        if (rName.empty()) {
            return r_empty_skin;
        }
        ModelPart& r_skin = mpModel->GetModelPart(rName);
        if (r_skin.NumberOfGeometries() == 0) {
            return r_skin;  // already a polyline of line conditions
        }
        if (!r_skin.HasSubModelPart("tessellation")) {
            TessellateNurbsSkin(r_skin, r_skin.CreateSubModelPart("tessellation"), tolerance);
        }
        return r_skin.GetSubModelPart("tessellation");
    };
    const ModelPart& r_outer_skin = get_polyline(outer_skin_name);
    const ModelPart& r_inner_skin = get_polyline(inner_skin_name);

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

    // Per skin (0: outer, 1: inner): its name and its identifier
    const std::array<std::string, 2> skin_names{outer_skin_name, inner_skin_name};
    const std::array<std::string, 2> skin_identifiers{"outer", "inner"};

    // First free node Id: node Ids are unique per root model part
    IndexType node_id = 1;
    for (const std::string& r_name : skin_names) {
        if (r_name.empty()) {
            continue;  // no skin of that kind, so nothing projects onto it
        }
        const ModelPart& r_root_model_part = mpModel->GetModelPart(r_name).GetRootModelPart();
        if (r_root_model_part.NumberOfNodes() > 0) {
            node_id = std::max(node_id, r_root_model_part.Nodes().back().Id() + 1);
        }
    }

    // Sub model part by name, created if missing
    const auto get_sub_model_part = [](ModelPart& rModelPart, const std::string& rName)
        -> ModelPart& {
        return rModelPart.HasSubModelPart(rName)
            ? rModelPart.GetSubModelPart(rName)
            : rModelPart.CreateSubModelPart(rName);
    };

    // Cut facets: points where alpha jumps, which get shifted boundary terms
    SizeType number_of_cut_facets = 0;

    for (auto& r_facet : r_facets.Conditions()) {
        const auto& r_minus = r_facet.GetGeometry().GetGeometryPart(0);
        const auto& r_plus = r_facet.GetGeometry().GetGeometryPart(1);

        if (r_minus.GetValue(VOLUME_FRACTION) == r_plus.GetValue(VOLUME_FRACTION)) {
            continue;  // same volume fraction on both sides: not cut, so no projection
        }
        ++number_of_cut_facets;

        // Closest point on each skin to the facet point, the same for both
        // sides, keeping the nearest (the outer one on ties)
        const array_1d<double, 3> point = r_minus.Center().Coordinates();
        IndexType skin = 0;
        NurbsProjectionResult closest = FindClosestSkinPoint(r_outer_skin, point);
        const NurbsProjectionResult inner = FindClosestSkinPoint(r_inner_skin, point);
        if (inner.SquaredDistance < closest.SquaredDistance) {
            closest = inner;
            skin = 1;
        }

        // <skin>.<layer>.projection_data for NURBS layers as SBM, else no layer
        ModelPart& r_skin = mpModel->GetModelPart(skin_names[skin]);
        ModelPart& r_projection_data = closest.LayerName.empty()
            ? get_sub_model_part(r_skin, "projection_data")
            : get_sub_model_part(
                get_sub_model_part(r_skin, closest.LayerName), "projection_data"
            );

        // Node at the closest skin point
        auto p_projection_node = r_projection_data.CreateNewNode(
            node_id++, closest.SkinCoordinates[0], closest.SkinCoordinates[1],
            closest.SkinCoordinates[2]
        );
        AssignBoundaryCurveGeometry(*p_projection_node, closest, skin == 0 ? 1.0 : -1.0);

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
        << number_of_cut_facets << " cut facets, on "
        << r_outer_skin.NumberOfConditions() + r_inner_skin.NumberOfConditions()
        << " skin segments." << std::endl;

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
        "skin_model_part_inner_initial_name" : "",
        "tessellation_tolerance"             : 1.0e-8
    })");
}

} // namespace Kratos
