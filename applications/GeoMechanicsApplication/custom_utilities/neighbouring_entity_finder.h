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

class KRATOS_API(GEO_MECHANICS_APPLICATION) NeighbouringEntityFinder
{
public:
    explicit NeighbouringEntityFinder(bool alsoSearchReverse = false);

    using BoundaryGeneratorByLocalDim = std::map<std::size_t, std::unique_ptr<BoundaryGenerator>>;

    template <typename EntityContainerType, typename CandidateNeighbourContainerType>
    void FindEntityNeighbours(EntityContainerType&             rEntities,
                              CandidateNeighbourContainerType& rCandidates,
                              BoundaryGeneratorByLocalDim&     rBoundaryGenerators)
    {
        for (int local_space_dimension = 0; local_space_dimension < 4; ++local_space_dimension) {
            NodeIdsToEntitiesHashMap map;
            if (!rBoundaryGenerators.contains(local_space_dimension)) continue;

            const auto& r_boundary_generator = *rBoundaryGenerators.at(local_space_dimension);
            for (auto& r_entity : rEntities) {
                if (r_entity.GetGeometry().LocalSpaceDimension() != local_space_dimension) continue;

                for (const auto& r_boundary_geometry : r_boundary_generator(r_entity.GetGeometry())) {
                    map.insert({GeometryUtilities::GetNodeIdsFromGeometry(r_boundary_geometry), {&r_entity}});
                }
            }

            if (map.empty()) continue;

            InitializeBoundaryMaps(map);
            FindEntityNeighboursBasedOnBoundaryType([&r_boundary_generator](const auto& rGeometry) {
                return r_boundary_generator(rGeometry);
            }, rCandidates);
        }
    }

private:
    void InitializeBoundaryMaps(NodeIdsToEntitiesHashMap GeometryNodeIdsToEntityMapping);
    void FindEntityNeighboursBasedOnBoundaryType(
        std::function<PointerVector<Geometry<Node>>(const Geometry<Node>&)> GenerateBoundaries,
        ModelPart::ElementsContainerType&                                   rElements);
    void AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
        Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries);
    void SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(const std::vector<std::size_t>& rConditionNodeIds,
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