//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/timer.h"
#include "utilities/entities_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "modeler/clean_up_problematic_triangles_modeler.h"

namespace Kratos
{

void CleanUpProblematicTrianglesModeler::SetupModelPart()
{
    KRATOS_TRY;

    Timer::Start("CleanUpProblematicTrianglesModeler");

    // We clean the mesh
    const std::string entity_type = mParameters["entity_type"].GetString();
    const IndexType first_node_id = mParameters["first_node_id"].GetInt();
    const IndexType first_element_id = mParameters["first_element_id"].GetInt();
    const IndexType first_condition_id = mParameters["first_condition_id"].GetInt();
    const double area_tolerance = mParameters["area_tolerance"].GetDouble();
    CleanUpProblematicGeometriesInMesh(*mpModelPart, entity_type, first_node_id, first_element_id, first_condition_id, area_tolerance);

    Timer::Stop("CleanUpProblematicTrianglesModeler");

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CleanUpProblematicTrianglesModeler::CleanUpProblematicGeometriesInMesh(
    ModelPart& rThisModelPart,
    const std::string& rEntityType,
    const IndexType FirstNodeId,
    const IndexType FirstElementId,
    const IndexType FirstConditionId,
    const double AreaTolerance
    )
{
    if (rEntityType == "element") {
        CleanUpProblematicGeometries<Element>(rThisModelPart, FirstNodeId, FirstElementId, AreaTolerance);
    } else if (rEntityType == "condition") {
        CleanUpProblematicGeometries<Condition>(rThisModelPart, FirstNodeId, FirstConditionId, AreaTolerance);
    } else if (rEntityType == "geometry") { // TODO: Implement geometry cleanup
        KRATOS_ERROR << "Geometry cleanup not yet supported. Sorry" << std::endl;
    } else {
        KRATOS_ERROR << "Unknown entity type: " << rEntityType << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <typename TEntityType>
void CleanUpProblematicTrianglesModeler::CleanUpProblematicGeometries(
    ModelPart& rThisModelPart,
    const IndexType FirstNodeId,
    const IndexType FirstEntityId,
    const double AreaTolerance
    )
{
    // Get the entities container
    auto& r_entities = EntitiesUtilities::GetEntities<TEntityType>(rThisModelPart);

    // Initialize variables
    double average_length = 0.0;
    double ref_area = 0.0;
    {
        // Compute the largest length of all geometries and the sum of lengths
        const double sum_length = block_for_each<SumReduction<double>>(r_entities, [](auto& rEntity) {
            const double length = rEntity.GetGeometry().Length();
            if (length > 0.0) {
                return length;
            } else {
                return 0.0;
            }
        });
        average_length = sum_length / static_cast<double>(r_entities.size());
        ref_area = average_length * average_length * AreaTolerance;
    }
    ref_area = ref_area * ref_area; // Using squared area tolerance to avoid expensive square roots, also will avoid numerical issues due to geometries that are 3 nodes in the same line, which squared area is negative and therefore area will be undefined (NaN)

    // Iterate until all null area triangles are removed
    std::size_t iter = 0;
    std::size_t null_area_triangles = ComputeNullAreaTriangles<TEntityType>(rThisModelPart, ref_area, AreaTolerance);
    KRATOS_INFO("CleanUpProblematicTrianglesModeler") << "Number of null area triangles: " << null_area_triangles << " in iteration " << iter << std::endl;
    while (null_area_triangles > 0) {
        // Initialize variables
        std::vector<typename TEntityType::Pointer> entities_to_remove;
        std::unordered_map<IndexType, std::unordered_set<IndexType>> replace_nodes;

        // Identify degenerated entities
        double distance_1, distance_2, distance_3;
        for (auto& r_entity : r_entities) {
            if (r_entity.IsNot(TO_ERASE)) {
                auto& r_geometry = r_entity.GetGeometry();

                // Check if the entity is a triangle
                if (r_geometry.GetGeometryType() != GeometryData::KratosGeometryType::Kratos_Triangle3D3) {
                    continue;
                }

                // Check that the area is below a certain tolerance
                const double squared_area = ComputeSquaredArea(r_geometry);
                if (squared_area < ref_area) {
                    // Get the nodes
                    const auto& r_node1 = r_geometry[0];
                    const auto& r_node2 = r_geometry[1];
                    const auto& r_node3 = r_geometry[2];

                    // Compute the squared distances
                    distance_1 = ComputeDistance(r_node1, r_node2);
                    distance_2 = ComputeDistance(r_node1, r_node3);
                    distance_3 = ComputeDistance(r_node2, r_node3);

                    // Compare which is the smallest of the three distances
                    const IndexType node_id_1 = r_node1.Id();
                    const IndexType node_id_2 = r_node2.Id();
                    const IndexType node_id_3 = r_node3.Id();
                    const double all_nodes_equal_tolerance = 1.0e-4;
                    if (distance_1/average_length < all_nodes_equal_tolerance && distance_2/average_length < all_nodes_equal_tolerance && distance_3/average_length < all_nodes_equal_tolerance) { // All nodes are the same
                        replace_nodes[node_id_1].insert(node_id_2);
                        replace_nodes[node_id_1].insert(node_id_3);
                        replace_nodes[node_id_2].insert(node_id_1);
                        replace_nodes[node_id_2].insert(node_id_3);
                        replace_nodes[node_id_3].insert(node_id_1);
                        replace_nodes[node_id_3].insert(node_id_2);
                    } else if (distance_1 < distance_2 && distance_1 < distance_3) {
                        // node1 and node2 are the same
                        replace_nodes[node_id_1].insert(node_id_2);
                        replace_nodes[node_id_2].insert(node_id_1);
                    } else if (distance_2 < distance_1 && distance_2 < distance_3) {
                        // node1 and node3 are the same
                        replace_nodes[node_id_1].insert(node_id_3);
                        replace_nodes[node_id_3].insert(node_id_1);
                    } else { // distance_3 is the smallest by elimination
                        // node2 and node3 are the same
                        replace_nodes[node_id_2].insert(node_id_3);
                        replace_nodes[node_id_3].insert(node_id_2);
                    }
                    r_entity.Set(TO_ERASE);
                }
            }
        }

        // Merge nodes to replace
        std::unordered_map<IndexType, IndexType> nodes_to_erase;
        for (auto& r_node_pair : replace_nodes) {
            const IndexType id = r_node_pair.first;
            auto& r_nodes_to_replace = r_node_pair.second;
            for (const IndexType node_id : r_nodes_to_replace) {
                const auto& r_sub_nodes_to_replace = replace_nodes[node_id];
                for (const IndexType sub_node_id : r_sub_nodes_to_replace) {
                    if (sub_node_id != id) {
                        r_nodes_to_replace.insert(sub_node_id);
                    }
                }
            }
            for (const IndexType node_id : r_nodes_to_replace) {
                replace_nodes.erase(node_id);
                nodes_to_erase.insert({node_id, id});
                rThisModelPart.GetNode(node_id).Set(TO_ERASE);
            }
        }

        // Clear replace nodes, not required anymore
        replace_nodes.clear();

        // Replace nodes in entities
        for (auto& r_entity : r_entities) {
            if (r_entity.IsNot(TO_ERASE)) {
                auto& r_geometry = r_entity.GetGeometry();
                // Iterate over the nodes of the geometry
                for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
                    const IndexType node_id = r_geometry[i].Id();
                    auto it_find_node = nodes_to_erase.find(node_id);
                    // If the node is in the list of nodes to replace
                    if (it_find_node != nodes_to_erase.end()) {
                        const IndexType replace_node_id = it_find_node->second;
                        // We check if the node is already on the geometry (will make null area triangles)
                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j) {
                            if (i != j && r_geometry[j].Id() == replace_node_id) {
                                r_entity.Set(TO_ERASE);
                            }
                        }
                        // Replace the node
                        if (r_entity.IsNot(TO_ERASE)) {
                            r_geometry(i) = rThisModelPart.pGetNode(replace_node_id);
                        }
                    }
                }
            }
        }

        // Recompute the number of null area triangles
        KRATOS_INFO("CleanUpProblematicTrianglesModeler") << "Number of null area triangles removed: " << null_area_triangles << " in iteration " << iter << std::endl;
        null_area_triangles = ComputeNullAreaTriangles<TEntityType>(rThisModelPart, ref_area, AreaTolerance);
        iter++;
        KRATOS_INFO_IF("CleanUpProblematicTrianglesModeler", null_area_triangles > 0) << "Number of null area triangles: " << null_area_triangles << " in iteration " << iter << std::endl;
    }

    // Remove entities and nodes marked to erase
    RemoveEntitiesAndNodes<TEntityType>(rThisModelPart);

    // Renumber
    {
        // Renumber nodes
        auto& r_nodes_array = rThisModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();
        IndexPartition<std::size_t>(r_nodes_array.size()).for_each([&](std::size_t i) {
            (it_node_begin + i)->SetId(i + FirstNodeId);
        });

        // Renumber entities
        const auto it_ent_begin = r_entities.begin();
        IndexPartition<std::size_t>(r_entities.size()).for_each([&](std::size_t i) {
            (it_ent_begin + i)->SetId(i + FirstEntityId);
        });
    }
}

// Explicit instantiation of template functions for Element and Condition
template void CleanUpProblematicTrianglesModeler::CleanUpProblematicGeometries<Element>(ModelPart& rThisModelPart, const IndexType FirstNodeId, const IndexType FirstEntityId, const double AreaTolerance);
template void CleanUpProblematicTrianglesModeler::CleanUpProblematicGeometries<Condition>(ModelPart& rThisModelPart, const IndexType FirstNodeId, const IndexType FirstEntityId, const double AreaTolerance);

/***********************************************************************************/
/***********************************************************************************/

// Helper function definitions

// RemoveEntitiesAndNodes specialization for Element
template <>
void CleanUpProblematicTrianglesModeler::RemoveEntitiesAndNodes<Element>(ModelPart& rThisModelPart)
{
    rThisModelPart.RemoveElements(TO_ERASE);
    rThisModelPart.RemoveNodes(TO_ERASE);
}

// RemoveEntitiesAndNodes specialization for Condition
template <>
void CleanUpProblematicTrianglesModeler::RemoveEntitiesAndNodes<Condition>(ModelPart& rThisModelPart)
{
    rThisModelPart.RemoveConditions(TO_ERASE);
    rThisModelPart.RemoveNodes(TO_ERASE);
}

/***********************************************************************************/
/***********************************************************************************/

double CleanUpProblematicTrianglesModeler::ComputeDistance(
    const Node& rNode1,
    const Node& rNode2
    )
{
    const array_1d<double, 3>& r_coords1 = rNode1.Coordinates();
    const array_1d<double, 3>& r_coords2 = rNode2.Coordinates();
    return norm_2(r_coords1 - r_coords2);
}

/***********************************************************************************/
/***********************************************************************************/

double CleanUpProblematicTrianglesModeler::ComputeSquaredArea(const GeometryType& rGeometry)
{
    const array_1d<double, 3>& r_point1 = rGeometry[0].Coordinates();
    const array_1d<double, 3>& r_point2 = rGeometry[1].Coordinates();
    const array_1d<double, 3>& r_point3 = rGeometry[2].Coordinates();

    const double a = MathUtils<double>::Norm3(r_point1 - r_point2);
    const double b = MathUtils<double>::Norm3(r_point2 - r_point3);
    const double c = MathUtils<double>::Norm3(r_point3 - r_point1);
    const double s = (a + b + c) / 2.0;
    return s * (s - a) * (s - b) * (s - c);
}

/***********************************************************************************/
/***********************************************************************************/

template <typename TEntityType>
std::size_t CleanUpProblematicTrianglesModeler::ComputeNullAreaTriangles(
    ModelPart& rThisModelPart,
    const double RefArea,
    const double AreaTolerance
    )
{
    const std::size_t null_area_triangles = block_for_each<SumReduction<std::size_t>>(EntitiesUtilities::GetEntities<TEntityType>(rThisModelPart), [&](auto& rEntity) {
        if (rEntity.IsNot(TO_ERASE)) {
            const auto& r_geometry = rEntity.GetGeometry();
            if (r_geometry.PointsNumber() == 3) {
                const double squared_area = ComputeSquaredArea(r_geometry);
                // Now check that the area is small enough
                if (squared_area < RefArea) {
                    return 1;
                }
            }
        }
        return 0;
    });
    return null_area_triangles;
}

// Explicit instantiation of template functions for Element and Condition
template std::size_t CleanUpProblematicTrianglesModeler::ComputeNullAreaTriangles<Element>(ModelPart& rThisModelPart, const double RefArea, const double AreaTolerance);
template std::size_t CleanUpProblematicTrianglesModeler::ComputeNullAreaTriangles<Condition>(ModelPart& rThisModelPart, const double RefArea, const double AreaTolerance);

}