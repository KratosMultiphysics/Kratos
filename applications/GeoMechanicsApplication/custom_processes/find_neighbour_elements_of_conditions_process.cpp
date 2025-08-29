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

namespace Kratos
{
void FindNeighbourElementsOfConditionsProcess::Execute()
{
    KRATOS_TRY

    // Next check that the conditions are oriented accordingly
    // to do so begin by putting all of the conditions in a map
    hashmap FacesMap;
    hashmap FacesMapSorted;

    for (auto itCond = mrModelPart.ConditionsBegin(); itCond != mrModelPart.ConditionsEnd(); ++itCond) {
        itCond->Set(VISITED, false);
        GeometryType& rGeometry = itCond->GetGeometry();

        DenseVector<IndexType> Ids(rGeometry.size());

        for (IndexType i = 0; i < Ids.size(); ++i) {
            rGeometry[i].Set(BOUNDARY, true);
            Ids[i] = rGeometry[i].Id();
        }

        // adds to the map
        FacesMap.insert(hashmap::value_type(Ids, std::vector<Condition::Pointer>({*itCond.base()})));

        DenseVector<int> IdsSorted = Ids;
        std::ranges::sort(IdsSorted);
        FacesMapSorted.insert(
            hashmap::value_type(IdsSorted, std::vector<Condition::Pointer>({*itCond.base()})));
    }

    if (FacesMap.empty()) return;

    // Now loop over all elements and check if one of the faces is in the "FacesMap"
    for (auto itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); ++itElem) {
        const auto& rGeometryElement    = itElem->GetGeometry();
        const auto  rBoundaryGeometries = rGeometryElement.GenerateBoundariesEntities();

        for (const auto& r_boundary_geometry : rBoundaryGeometries) {
            DenseVector<IndexType> FaceIds(r_boundary_geometry.size());
            std::ranges::transform(r_boundary_geometry, FaceIds.begin(),
                                   [](const Node& rNode) { return rNode.Id(); });

            hashmap::iterator itFace = FacesMap.find(FaceIds);
            if (itFace == FacesMap.end() && (rGeometryElement.LocalSpaceDimension() == 3 ||
                                             (rGeometryElement.LocalSpaceDimension() == 2 &&
                                              dynamic_cast<const InterfaceElement*>(&(*itElem))))) {
                // condition is not found but might be a problem of ordering in 3D geometries!
                DenseVector<int> FaceIdsSorted = FaceIds;
                std::sort(FaceIdsSorted.begin(), FaceIdsSorted.end());
                hashmap::iterator itFaceSorted = FacesMapSorted.find(FaceIdsSorted);
                if (itFaceSorted != FacesMapSorted.end()) {
                    // try different orderings
                    using enum GeometryData::KratosGeometryOrderType;
                    if (r_boundary_geometry.GetGeometryOrderType() == Kratos_Linear_Order)
                        itFace = FindPermutations(FaceIds, FacesMap);
                    else if (r_boundary_geometry.GetGeometryOrderType() == Kratos_Quadratic_Order)
                        itFace = FindPermutationsQuadratic(FaceIds, FacesMap);
                }
            }

            if (itFace != FacesMap.end()) {
                // condition is found!
                // but check if there are more than one condition on the element
                CheckForMultipleConditionsOnElement(FacesMap, itFace, itElem);
            }
        }
    }

    // check that all of the conditions belong to at least an element.
    bool all_conditions_visited = CheckIfAllConditionsAreVisited();

    if (all_conditions_visited) {
        // if all conditions are found, no need for further checks:
        return;
    } else {
        // Now try point loads:
        for (auto itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); ++itElem) {
            const auto& rGeometryElement = itElem->GetGeometry();
            const auto  rPointGeometries = rGeometryElement.GeneratePoints();

            for (IndexType iPoint = 0; iPoint < rPointGeometries.size(); ++iPoint) {
                DenseVector<IndexType> PointIds(rPointGeometries[iPoint].size());

                // Points
                for (IndexType iNode = 0; iNode < PointIds.size(); ++iNode) {
                    PointIds[iNode] = rPointGeometries[iPoint][iNode].Id();
                }

                hashmap::iterator itFace = FacesMap.find(PointIds);
                if (itFace != FacesMap.end()) {
                    // condition is found!
                    // but check if there are more than one condition on the element
                    CheckForMultipleConditionsOnElement(FacesMap, itFace, itElem);
                }
            }
        }
    }

    // check that all of the conditions belong to at least an element.
    all_conditions_visited = CheckIfAllConditionsAreVisited();

    if (all_conditions_visited) {
        // if all conditions are found, no need for further checks:
        return;
    } else {
        // check edges of 3D geometries:
        // Now loop over all elements and check if one of the faces is in the "FacesMap"
        for (auto itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); ++itElem) {
            const auto& rGeometryElement = itElem->GetGeometry();
            if (rGeometryElement.LocalSpaceDimension() == 3) {
                const auto rBoundaryGeometries = rGeometryElement.GenerateEdges();

                for (IndexType iEdge = 0; iEdge < rBoundaryGeometries.size(); ++iEdge) {
                    DenseVector<IndexType> EdgeIds(rBoundaryGeometries[iEdge].size());

                    // edges for 3D elements
                    for (IndexType iNode = 0; iNode < EdgeIds.size(); ++iNode) {
                        EdgeIds[iNode] = rBoundaryGeometries[iEdge][iNode].Id();
                    }

                    hashmap::iterator itFace = FacesMap.find(EdgeIds);
                    // There might be a need to check this for different types of 3D elements
                    // as the ordering numbers might be inconsistent

                    if (itFace != FacesMap.end()) {
                        // condition is found!
                        // but check if there are more than one condition on the element
                        CheckForMultipleConditionsOnElement(FacesMap, itFace, itElem);
                    }
                }
            }
        }
    }

    // check that all of the conditions belong to at least an element.
    all_conditions_visited = CheckIfAllConditionsAreVisited();

    if (all_conditions_visited) {
        // if all conditions are found, no need for further checks:
        return;
    }

    // check 1D elements, note that this has to happen after procedures to find 2 and 3d neighbours are alredy performed, such that 1D elements are only added
    // as neighbours when the condition is not neighbouring 2D or 3D elements
    this->CheckIf1DElementIsNeighbour(FacesMap);

    // check that all of the conditions belong to at least an element. Throw an error otherwise (this is particularly useful in mpi)
    all_conditions_visited = true;
    for (const auto& rCond : mrModelPart.Conditions()) {
        if (rCond.IsNot(VISITED)) {
            all_conditions_visited = false;
            KRATOS_INFO("Condition without any corresponding element, ID ") << rCond.Id() << std::endl;
        }
    }
    KRATOS_ERROR_IF_NOT(all_conditions_visited)
        << "Some conditions found without any corresponding element" << std::endl;

    KRATOS_CATCH("")
}

bool FindNeighbourElementsOfConditionsProcess::CheckIfAllConditionsAreVisited() const
{
    return std::ranges::all_of(mrModelPart.Conditions(),
                               [](const auto& r_cond) { return r_cond.Is(VISITED); });
}

void FindNeighbourElementsOfConditionsProcess::CheckIf1DElementIsNeighbour(hashmap& rFacesMap)
{
    // Now loop over all elements and check if one of the faces is in the "FacesMap"
    for (auto itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); ++itElem) {
        const auto& r_geometry_element = itElem->GetGeometry();

        // for 1D elements, the edge geometry is the same as the element geometry
        if (r_geometry_element.LocalSpaceDimension() == 1) {
            const auto rBoundaryGeometries = PointerVector(r_geometry_element.GenerateEdges());

            for (IndexType iFace = 0; iFace < rBoundaryGeometries.size(); ++iFace) {
                DenseVector<int> FaceIds(rBoundaryGeometries[iFace].size());

                const auto& r_nodes = rBoundaryGeometries[iFace];

                // get face node IDs
                std::transform(r_nodes.begin(), r_nodes.end(), FaceIds.begin(),
                               [](const auto& r_node) { return r_node.Id(); });

                hashmap::iterator itFace = rFacesMap.find(FaceIds);

                if (itFace != rFacesMap.end()) {
                    // condition is found!
                    // but check if there are more than one condition on the element
                    CheckForMultipleConditionsOnElement(rFacesMap, itFace, itElem);
                }
            }
        }
    }
}

void FindNeighbourElementsOfConditionsProcess::CheckForMultipleConditionsOnElement(
    hashmap& rFacesMap, hashmap::iterator& rItFace, PointerVector<Element>::iterator pItElem)
{
    const std::pair<hashmap::iterator, hashmap::iterator> face_pair = rFacesMap.equal_range(rItFace->first);
    for (hashmap::iterator it = face_pair.first; it != face_pair.second; ++it) {
        std::vector<Condition::Pointer>& r_conditions = it->second;

        GlobalPointersVector<Element> vector_of_neighbours;
        vector_of_neighbours.resize(1);
        vector_of_neighbours(0) = Element::WeakPointer(*pItElem.base());

        for (Condition::Pointer p_condition : r_conditions) {
            p_condition->Set(VISITED, true);
            p_condition->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
        }
    }
}

hashmap::iterator FindNeighbourElementsOfConditionsProcess::FindPermutations(DenseVector<int> FaceIds,
                                                                             hashmap& FacesMap) const
{
    for (std::size_t i = 0; i < FaceIds.size() - 1; ++i) {
        std::ranges::rotate(FaceIds, FaceIds.begin() + 1);

        hashmap::iterator itFace = FacesMap.find(FaceIds);
        if (itFace != FacesMap.end()) return itFace;
    }
    return FacesMap.end();
}

hashmap::iterator FindNeighbourElementsOfConditionsProcess::FindPermutationsQuadratic(DenseVector<int> FaceIds,
                                                                                      hashmap& FacesMap) const
{
    for (std::size_t i = 0; i < FaceIds.size() / 2 - 1; ++i) {
        std::rotate(FaceIds.begin(), FaceIds.begin() + 1, FaceIds.begin() + FaceIds.size() / 2);
        std::rotate(FaceIds.begin() + FaceIds.size() / 2, FaceIds.begin() + FaceIds.size() / 2 + 1,
                    FaceIds.end());

        hashmap::iterator itFace = FacesMap.find(FaceIds);
        if (itFace != FacesMap.end()) return itFace;
    }
    return FacesMap.end();
}
} // namespace Kratos