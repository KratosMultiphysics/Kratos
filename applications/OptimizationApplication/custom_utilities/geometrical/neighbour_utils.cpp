//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <vector>
#include <unordered_map>

// Project includes
#include "containers/global_pointers_vector.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "neighbour_utils.h"

namespace Kratos {

void NeighbourUtils::InitializeParentElementsForConditions(ModelPart& rModelPart)
{
    KRATOS_TRY

    using key_type = std::vector<IndexType>;
    using map_type = std::unordered_map<key_type, Condition*, KeyHasherRange<key_type>, KeyComparorRange<key_type>>;

    KRATOS_ERROR_IF(rModelPart.NumberOfElements() == 0)
        << rModelPart.FullName() << " does not contain any elements which are required for initialize parent elements for conditions.";

    KRATOS_ERROR_IF(rModelPart.NumberOfConditions() == 0)
        << rModelPart.FullName() << " does not contain any conditions which are required for initialize parent elements for conditions.";

    // first generate the map of condition and its node ids
    const auto& map = block_for_each<MapReduction<map_type>>(rModelPart.Conditions(), [&](auto& rCondition) {
        key_type node_ids;
        for (const auto& r_node : rCondition.GetGeometry()) {
            node_ids.push_back(r_node.Id());
        }
        std::sort(node_ids.begin(), node_ids.end());

        // we reset the NEIGHBOUR_ELEMENTS vector of the corresponding condition
        GlobalPointersVector<Element> vector_of_neighbours;
        rCondition.SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);

        return std::make_pair<key_type, Condition*>(std::forward<key_type>(node_ids), &rCondition);
    });

    // now iterate through element faces, edges, nodes
    block_for_each(rModelPart.Elements(), [&](auto& rElement) {
        const auto& r_geometry = rElement.GetGeometry();
        const IndexType dimension = r_geometry.LocalSpaceDimension();

        if (dimension == 3) {
            // 3 dimensional geometries can have conditions which can be either face, edge or point.
            // hence faces, edges and points needs to be checked.

            // check for faces.
            const auto& r_boundary_faces = r_geometry.GenerateFaces();
            for (const auto& r_boundary_face : r_boundary_faces) {
                key_type node_ids;
                for (const auto& r_boundary_node : r_boundary_face) {
                    node_ids.push_back(r_boundary_node.Id());
                }
                std::sort(node_ids.begin(), node_ids.end());

                const auto& p_found_it = map.find(node_ids);
                if (p_found_it != map.end()) {
                    // this is done without any critical section because, for each element, there can be only
                    // one face shared among condition and element. So there won't be any race conditions.
                    auto& r_neighbour_elements_vector = p_found_it->second->GetValue(NEIGHBOUR_ELEMENTS);
                    r_neighbour_elements_vector.push_back(Element::WeakPointer(&rElement));
                }
            }

            // check for edges.
            const auto& r_boundary_edges = r_geometry.GenerateEdges();
            for (const auto& r_boundary_edge : r_boundary_edges) {
                key_type node_ids;
                for (const auto& r_boundary_node : r_boundary_edge) {
                    node_ids.push_back(r_boundary_node.Id());
                }
                std::sort(node_ids.begin(), node_ids.end());

                const auto& p_found_it = map.find(node_ids);
                if (p_found_it != map.end()) {
                    KRATOS_CRITICAL_SECTION
                    auto& r_neighbour_elements_vector = p_found_it->second->GetValue(NEIGHBOUR_ELEMENTS);
                    r_neighbour_elements_vector.push_back(Element::WeakPointer(&rElement));
                }
            }
        } else if (dimension == 2) {
            // 2 dimensional geometries can have conditions which can be either edge or point.
            // hence edges and points needs to be checked.

            // check for edges.
            const auto& r_boundary_edges = r_geometry.GenerateEdges();
            for (const auto& r_boundary_edge : r_boundary_edges) {
                key_type node_ids;
                for (const auto& r_boundary_node : r_boundary_edge) {
                    node_ids.push_back(r_boundary_node.Id());
                }
                std::sort(node_ids.begin(), node_ids.end());

                const auto& p_found_it = map.find(node_ids);
                if (p_found_it != map.end()) {
                    // this is done without any critical section because, for each element, there can be only
                    // one edge shared among condition and element. So there won't be any race conditions.
                    auto& r_neighbour_elements_vector = p_found_it->second->GetValue(NEIGHBOUR_ELEMENTS);
                    r_neighbour_elements_vector.push_back(Element::WeakPointer(&rElement));
                }
            }
        } else {
            KRATOS_ERROR << "Unsupported LocalSpaceDimension for element with id "
                         << rElement.Id() << " found in " << rModelPart.FullName()
                         << ". [ dimension = " << dimension << " ].\n";
        }

        // check for points.
        const auto& r_boundary_points = r_geometry.GeneratePoints();
        for (const auto& r_boundary_point : r_boundary_points) {
            const auto& p_found_it = map.find({r_boundary_point[0].Id()});
            if (p_found_it != map.end()) {
                KRATOS_CRITICAL_SECTION
                auto& r_neighbour_elements_vector = p_found_it->second->GetValue(NEIGHBOUR_ELEMENTS);
                r_neighbour_elements_vector.push_back(Element::WeakPointer(&rElement));
            }
        }
    });

    KRATOS_CATCH("");
}

std::unordered_map<IndexType, std::vector<IndexType>> NeighbourUtils::GetConditionIdAndParentElementIdsMap(const ModelPart& rModelPart)
{
    KRATOS_TRY

    return block_for_each<MapReduction<std::unordered_map<IndexType, std::vector<IndexType>>>>(rModelPart.Conditions(), [&](const auto& rCondition) {
        std::vector<IndexType> ids;
        for (const auto& p_element : rCondition.GetValue(NEIGHBOUR_ELEMENTS).GetContainer()) {
            ids.push_back(p_element->Id());
        }
        std::sort(ids.begin(), ids.end());
        return std::make_pair(rCondition.Id(), ids);
    });

    KRATOS_CATCH("");
}

} // namespace Kratos