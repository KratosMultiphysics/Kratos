// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "boundary_generators.h"
#include "geometry_utilities.h"
#include "includes/geometrical_object.h"
#include "includes/key_hash.h"
#include "includes/kratos_export_api.h"
#include "includes/model_part.h"

#include <unordered_map>

namespace Kratos
{

using NodeIdsToEntitiesHashMap       = std::unordered_multimap<std::vector<std::size_t>,
                                                               std::vector<GeometricalObject::Pointer>,
                                                               KeyHasherRange<std::vector<std::size_t>>,
                                                               KeyComparorRange<std::vector<std::size_t>>>;
using SortedToUnsortedNodeIdsHashMap = std::unordered_multimap<std::vector<std::size_t>,
                                                               std::vector<std::size_t>,
                                                               KeyHasherRange<std::vector<std::size_t>>,
                                                               KeyComparorRange<std::vector<std::size_t>>>;

/**
 * @class NeighbouringElementFinder
 * @brief Utility class for finding neighbouring elements of generic geometrical objects based on shared boundary geometries.
 * @details This class identifies neighbour relationships between geometrical entities by comparing their boundary
 * geometries (edges, faces, etc.). It supports both forward and reverse search modes, and can handle linear and
 * quadratic elements with potentially rotated or permuted node orderings.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) NeighbouringElementFinder
{
public:
    /**
     * @brief Constructor for NeighbouringEntityFinder
     * @param alsoSearchReverse If true, also searches for reversed node orderings when matching boundaries.
     * This is useful when elements may have opposite orientations but share the same boundary.
     */
    explicit NeighbouringElementFinder(bool alsoSearchReverse = false);

    using BoundaryGeneratorByLocalDim = std::map<std::size_t, std::unique_ptr<BoundaryGenerator>>;

    /**
     * @brief Finds and assigns neighbour elements based on shared boundary geometries.
     * @details This method iterates through all local space dimensions (0-3) and uses the appropriate boundary
     * generators to extract boundary geometries (e.g., edges for 2D elements, faces for 3D elements). It then
     * compares these boundaries with those of candidate elements to identify neighbours. Matched elements are
     * stored in the NEIGHBOUR_ELEMENTS variable of the entity.
     *
     * The method supports both linear and quadratic elements and can handle node orderings that are rotated or
     * permuted, as long as the node IDs of shared boundaries match.
     *
     * @tparam EntityContainerType Container type for entities (e.g., ModelPart::ElementsContainerType,
     * ModelPart::ConditionsContainerType)
     * @param rEntities Container of entities for which to find neighbours
     * @param rCandidateElements Container of candidate elements that may be neighbours
     * @param rBoundaryGenerators Map of local space dimensions to boundary generators (e.g., EdgesGenerator for
     * 2D, FacesGenerator for 3D). The key is the local space dimension and the value is the generator.
     */
    template <typename EntityContainerType>
    void FindEntityNeighbours(EntityContainerType&              rEntities,
                              ModelPart::ElementsContainerType& rCandidateElements,
                              BoundaryGeneratorByLocalDim&      rBoundaryGenerators)
    {
        for (std::size_t local_space_dimension = 0; local_space_dimension < 4; ++local_space_dimension) {
            if (!rBoundaryGenerators.contains(local_space_dimension)) continue;

            NodeIdsToEntitiesHashMap map;
            const auto& r_boundary_generator = *rBoundaryGenerators.at(local_space_dimension);
            for (auto& r_entity : rEntities) {
                if (r_entity.GetGeometry().LocalSpaceDimension() != local_space_dimension) continue;

                for (const auto& r_boundary_geometry : r_boundary_generator(r_entity.GetGeometry())) {
                    map.insert({GeometryUtilities::GetNodeIdsFromGeometry(r_boundary_geometry), {&r_entity}});
                }
            }

            if (map.empty()) continue;

            InitializeBoundaryMaps(map);
            FindEntityNeighboursBasedOnBoundaryType(r_boundary_generator, rCandidateElements);
        }
    }

private:
    void InitializeBoundaryMaps(NodeIdsToEntitiesHashMap GeometryNodeIdsToEntityMapping);
    void FindEntityNeighboursBasedOnBoundaryType(const BoundaryGenerator& rBoundaryGenerator,
                                                 ModelPart::ElementsContainerType& rElements);
    void AddNeighbouringElementsToEntitiesBasedOnOverlappingBoundaryGeometries(
        Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries);
    void SetElementAsNeighbourOfAllEntitiesWithIdenticalNodeIds(const std::vector<std::size_t>& rEntityNodeIds,
                                                                Element* pElement);
    void        SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(Element& rElement,
                                                                   const std::vector<std::size_t>& rElementBoundaryNodeIds,
                                                                   GeometryData::KratosGeometryOrderType OrderType);
    static bool AreRotatedEquivalents(const std::vector<std::size_t>&       rFirst,
                                      const std::vector<std::size_t>&       rSecond,
                                      GeometryData::KratosGeometryOrderType OrderType);
    static bool AreLinearRotatedEquivalents(std::vector<std::size_t>        First,
                                            const std::vector<std::size_t>& rSecond);
    static bool AreQuadraticRotatedEquivalents(std::vector<std::size_t>        First,
                                               const std::vector<std::size_t>& rSecond);

    NodeIdsToEntitiesHashMap       mGeometryNodeIdsToEntities;
    SortedToUnsortedNodeIdsHashMap mSortedToUnsortedEntityNodeIds;
    bool                           mAlsoSearchReverse = false;
};
} // namespace Kratos