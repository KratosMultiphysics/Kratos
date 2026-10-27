//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License (kratos/license.txt)

// System includes
#include <algorithm>
#include <cmath>

// Project includes
#include "utilities/multipatch_utilities.h"

namespace Kratos
{

FemMultipatchUtilities::FemMultipatchUtilities(ModelPart& rBackground)
    : mrBackground(rBackground)
{
    RebuildTopology();
}

void FemMultipatchUtilities::RebuildTopology()
{
    KRATOS_TRY

    mElements.clear();
    mEdges.clear();
    mElements.reserve(mrBackground.NumberOfElements());
    std::map<NodePairType, std::size_t> edge_indices;
    for (auto& r_element : mrBackground.Elements()) {
        auto& r_geometry = r_element.GetGeometry();
        const auto type = r_geometry.GetGeometryType();
        KRATOS_ERROR_IF(type != GeometryData::KratosGeometryType::Kratos_Triangle2D3 &&
                        type != GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
            << "FEM multipatch background requires P1 triangles or Q1 quadrilaterals."
            << std::endl;
        auto p_element = mrBackground.pGetElement(r_element.Id());
        mElements.push_back(p_element);
        // Use geometry topology; quadrilateral diagonals are never edges.
        for (const auto& r_edge : r_geometry.GenerateEdges()) {
            const NodePairType nodes{r_edge[0].Id(), r_edge[1].Id()};
            auto key = nodes;
            std::sort(key.begin(), key.end());
            const auto [it, inserted] = edge_indices.emplace(key, mEdges.size());
            if (inserted) {
                mEdges.emplace_back();
            }
            auto& r_data = mEdges[it->second];
            const std::size_t side = r_data.Elements[0] ? 1 : 0;
            KRATOS_ERROR_IF(r_data.Elements[side])
                << "Non-manifold background edge " << key[0] << ", " << key[1] << std::endl;
            r_data.Elements[side] = p_element;
            r_data.Nodes[side] = nodes;
        }
    }

    KRATOS_CATCH("")
}

FemMultipatchUtilities::OwnerMapType FemMultipatchUtilities::UpdateBackgroundActivation(
    const std::vector<BoxType>& rBoxes,
    const double Tolerance)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(!std::isfinite(Tolerance) || Tolerance <= 0.0)
        << "Tolerance must be finite and positive." << std::endl;
    OwnerMapType inactive_owners;
    for (auto& p_element : mElements) {
        p_element->Set(ACTIVE, true);
        for (std::size_t i = 0; i < rBoxes.size(); ++i) {
            const auto& r_box = rBoxes[i];
            // All vertices must lie in the box. Centroid classification would
            // silently accept cut elements, which Phase 1 does not support.
            const auto& r_geometry = p_element->GetGeometry();
            const bool inside = std::all_of(r_geometry.begin(), r_geometry.end(),
                [&r_box, Tolerance](const Node& rNode) {
                    return rNode.X() >= r_box[0] - Tolerance &&
                           rNode.Y() >= r_box[1] - Tolerance &&
                           rNode.X() <= r_box[2] + Tolerance &&
                           rNode.Y() <= r_box[3] + Tolerance;
                });
            if (inside) {
                p_element->Set(ACTIVE, false);
                inactive_owners.emplace(p_element->Id(), i + 1);
                break;
            }
        }
    }
    return inactive_owners;

    KRATOS_CATCH("")
}

FemMultipatchUtilities::SurrogateEdgesType FemMultipatchUtilities::ExtractSurrogateEdges(
    const OwnerMapType& rInactiveOwners) const
{
    KRATOS_TRY

    SurrogateEdgesType surrogate_edges;
    for (const auto& r_edge : mEdges) {
        if (!r_edge.Elements[1]) {
            continue; // Exclude the outer background boundary.
        }
        const bool first_active = r_edge.Elements[0]->Is(ACTIVE);
        const bool second_active = r_edge.Elements[1]->Is(ACTIVE);
        if (first_active != second_active) {
            const std::size_t active_side = first_active ? 0 : 1;
            const auto inactive_id = r_edge.Elements[1 - active_side]->Id();
            const auto it = rInactiveOwners.find(inactive_id);
            KRATOS_ERROR_IF(it == rInactiveOwners.end())
                << "Missing internal patch owner for inactive element " << inactive_id << std::endl;
            // Preserve the active element's own edge orientation.
            surrogate_edges[it->second].push_back(r_edge.Nodes[active_side]);
        }
    }
    return surrogate_edges;

    KRATOS_CATCH("")
}

FemMultipatchUtilities::GaussPointsType FemMultipatchUtilities::GetGaussPoints(
    const ModelPart& rModelPart,
    const bool ActiveOnly)
{
    KRATOS_TRY

    GaussPointsType result;
    for (const auto& r_element : rModelPart.Elements()) {
        if (ActiveOnly && !r_element.Is(ACTIVE)) {
            continue;
        }
        const auto& r_geometry = r_element.GetGeometry();
        const auto& r_points = r_geometry.IntegrationPoints(r_element.GetIntegrationMethod());
        for (const auto& r_point : r_points) {
            array_1d<double, 3> global_coordinates = ZeroVector(3);
            r_geometry.GlobalCoordinates(global_coordinates, r_point);
            result.push_back({r_element.Id(),
                {global_coordinates[0], global_coordinates[1], global_coordinates[2]}});
        }
    }
    return result;

    KRATOS_CATCH("")
}

std::vector<std::tuple<FemMultipatchUtilities::IndexType, std::array<double, 3>, double>>
FemMultipatchUtilities::GetGaussPointsWithWeights(const ModelPart& rModelPart, const bool ActiveOnly)
{
    std::vector<std::tuple<IndexType, std::array<double, 3>, double>> result;
    for (const auto& e : rModelPart.Elements()) {
        if (ActiveOnly && e.IsDefined(ACTIVE) && e.IsNot(ACTIVE)) continue;
        const auto& geometry=e.GetGeometry();
        const auto method=e.GetIntegrationMethod();
        const auto& points=geometry.IntegrationPoints(method);
        for (std::size_t i=0;i<points.size();++i) {
            Geometry<Node>::CoordinatesArrayType global=ZeroVector(3);
            geometry.GlobalCoordinates(global,points[i]);
            result.emplace_back(e.Id(),std::array<double,3>{global[0],global[1],global[2]},
                points[i].Weight()*geometry.DeterminantOfJacobian(i,method));
        }
    }
    return result;
}

} // namespace Kratos
