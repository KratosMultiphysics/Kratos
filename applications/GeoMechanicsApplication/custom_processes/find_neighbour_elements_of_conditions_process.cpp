// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#include "custom_processes/find_neighbour_elements_of_conditions_process.hpp"

#include "custom_elements/interface_element.h"
#include "geometries/geometry.h"
#include "includes/kratos_flags.h"
#include <memory>

namespace Kratos
{

void FindNeighbourElementsOfConditionsProcess::Execute()
{
    KRATOS_TRY

    if (mrModelPart.Conditions().empty()) return;

    hashmap condition_node_ids_to_condition;
    for (auto& r_condition : mrModelPart.Conditions()) {
        r_condition.Set(VISITED, false);
        auto& r_geometry = r_condition.GetGeometry();

        std::vector<IndexType> Ids(r_geometry.size());
        std::ranges::transform(r_geometry, Ids.begin(), [](const auto& rNode) { return rNode.Id(); });
        for (auto& rNode : r_geometry)
            rNode.Set(BOUNDARY, true);

        condition_node_ids_to_condition.insert(hashmap::value_type(Ids, {&r_condition}));
    }

    hashmap sorted_condition_node_ids_to_condition;
    std::ranges::transform(condition_node_ids_to_condition, std::inserter(sorted_condition_node_ids_to_condition, sorted_condition_node_ids_to_condition.end()), [](const auto& rPair) {
        auto sorted_ids = rPair.first;
        std::ranges::sort(sorted_ids);
        return std::make_pair(sorted_ids, rPair.second);
    });

    for (auto& r_element : mrModelPart.Elements()) {
        const auto& rGeometryElement    = r_element.GetGeometry();
        const auto  rBoundaryGeometries = rGeometryElement.GenerateBoundariesEntities();
        AddNeighboringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
            condition_node_ids_to_condition, sorted_condition_node_ids_to_condition, r_element, rBoundaryGeometries);
    }

    if (AllConditionsAreVisited()) return;

    // Now try point loads:
    for (auto& r_element : mrModelPart.Elements()) {
        const auto& rGeometryElement    = r_element.GetGeometry();
        const auto  rBoundaryGeometries = rGeometryElement.GeneratePoints();

        AddNeighboringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
            condition_node_ids_to_condition, sorted_condition_node_ids_to_condition, r_element, rBoundaryGeometries);
    }

    if (AllConditionsAreVisited()) return;

    // check edges of 3D geometries:
    // Now loop over all elements and check if one of the faces is in the "FacesMap"
    for (auto& r_element : mrModelPart.Elements()) {
        const auto& r_geometry_element = r_element.GetGeometry();
        if (r_geometry_element.LocalSpaceDimension() == 3) {
            const auto& r_boundary_geometries = r_geometry_element.GenerateEdges();

            AddNeighboringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
                condition_node_ids_to_condition, sorted_condition_node_ids_to_condition, r_element, r_boundary_geometries);
        }
    }

    if (AllConditionsAreVisited()) return;

    // check 1D elements, note that this has to happen after procedures to find 2 and 3d neighbours are already performed, such that 1D elements are only added
    // as neighbours when the condition is not neighbouring 2D or 3D elements
    this->CheckIf1DElementIsNeighbour(condition_node_ids_to_condition);

    // check that all of the conditions belong to at least an element. Throw an error otherwise (this is particularly useful in mpi)
    auto all_conditions_visited = true;
    for (const auto& r_condition : mrModelPart.Conditions()) {
        if (r_condition.IsNot(VISITED)) {
            all_conditions_visited = false;
            KRATOS_INFO("Condition without any corresponding element, ID ") << r_condition.Id() << std::endl;
        }
    }
    KRATOS_ERROR_IF_NOT(all_conditions_visited)
        << "Some conditions found without any corresponding element" << std::endl;

    KRATOS_CATCH("")
}

void FindNeighbourElementsOfConditionsProcess::AddNeighboringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
    hashmap& FacesMap, const hashmap& FacesMapSorted, Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries)
{
    for (const auto& r_boundary_geometry : rBoundaryGeometries) {
        std::vector<IndexType> face_ids(r_boundary_geometry.size());
        std::ranges::transform(r_boundary_geometry, face_ids.begin(),
                               [](const Node& rNode) { return rNode.Id(); });

        auto itFace = FacesMap.find(face_ids);
        if (itFace == FacesMap.end() && r_boundary_geometry.LocalSpaceDimension() == 2) {
            // condition is not found but might be a problem of ordering in 2D boundary geometries!
            std::vector<std::size_t> face_ids_sorted = face_ids;
            std::ranges::sort(face_ids_sorted);
            if (FacesMapSorted.contains(face_ids_sorted)) {
                switch (r_boundary_geometry.GetGeometryOrderType()) {
                    using enum GeometryData::KratosGeometryOrderType;
                case Kratos_Linear_Order:
                    itFace = FindPermutations(face_ids, FacesMap);
                    break;
                case Kratos_Quadratic_Order:
                    itFace = FindPermutationsQuadratic(face_ids, FacesMap);
                    break;
                default:
                    break;
                }
            }
        }

        if (itFace != FacesMap.end()) {
            // condition is found!
            // but check if there are more than one condition on the element
            CheckForMultipleConditionsOnElement(FacesMap, itFace, &rElement);
        }
    }
}

bool FindNeighbourElementsOfConditionsProcess::AllConditionsAreVisited() const
{
    return std::ranges::all_of(mrModelPart.Conditions(),
                               [](const auto& rCondition) { return rCondition.Is(VISITED); });
}

void FindNeighbourElementsOfConditionsProcess::CheckIf1DElementIsNeighbour(hashmap& rFacesMap)
{
    // Now loop over all elements and check if one of the faces is in the "FacesMap"
    for (auto& r_element : mrModelPart.Elements()) {
        const auto& r_geometry_element = r_element.GetGeometry();

        // for 1D elements, the edge geometry is the same as the element geometry
        if (r_geometry_element.LocalSpaceDimension() == 1) {
            const auto boundary_geometries = PointerVector(r_geometry_element.GenerateEdges());

            for (IndexType iFace = 0; iFace < boundary_geometries.size(); ++iFace) {
                std::vector<std::size_t> FaceIds(boundary_geometries[iFace].size());

                const auto& r_nodes = boundary_geometries[iFace];

                // get face node IDs
                std::transform(r_nodes.begin(), r_nodes.end(), FaceIds.begin(),
                               [](const auto& r_node) { return r_node.Id(); });

                auto itFace = rFacesMap.find(FaceIds);

                if (itFace != rFacesMap.end()) {
                    // condition is found!
                    // but check if there are more than one condition on the element
                    CheckForMultipleConditionsOnElement(rFacesMap, itFace, &r_element);
                }
            }
        }
    }
}

void FindNeighbourElementsOfConditionsProcess::CheckForMultipleConditionsOnElement(hashmap& rFacesMap,
                                                                                   const hashmap::iterator& rItFace,
                                                                                   Element* pElement)
{
    const auto face_pair = rFacesMap.equal_range(rItFace->first);
    for (auto it = face_pair.first; it != face_pair.second; ++it) {
        auto&                         r_conditions = it->second;
        GlobalPointersVector<Element> vector_of_neighbours;
        vector_of_neighbours.resize(1);
        vector_of_neighbours(0) = Element::WeakPointer(pElement);

        for (auto& p_condition : r_conditions) {
            p_condition->Set(VISITED, true);
            p_condition->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
        }
    }
}

hashmap::iterator FindNeighbourElementsOfConditionsProcess::FindPermutations(std::vector<std::size_t> FaceIds,
                                                                             hashmap& FacesMap) const
{
    for (std::size_t i = 0; i < FaceIds.size() - 1; ++i) {
        std::ranges::rotate(FaceIds, FaceIds.begin() + 1);

        auto itFace = FacesMap.find(FaceIds);
        if (itFace != FacesMap.end()) return itFace;
    }
    return FacesMap.end();
}

hashmap::iterator FindNeighbourElementsOfConditionsProcess::FindPermutationsQuadratic(std::vector<std::size_t> FaceIds,
                                                                                      hashmap& FacesMap) const
{
    for (std::size_t i = 0; i < FaceIds.size() / 2 - 1; ++i) {
        std::rotate(FaceIds.begin(), FaceIds.begin() + 1, FaceIds.begin() + FaceIds.size() / 2);
        std::rotate(FaceIds.begin() + FaceIds.size() / 2, FaceIds.begin() + FaceIds.size() / 2 + 1,
                    FaceIds.end());

        auto itFace = FacesMap.find(FaceIds);
        if (itFace != FacesMap.end()) return itFace;
    }
    return FacesMap.end();
}
} // namespace Kratos