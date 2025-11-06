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
//                   Richard Faasse
//

#pragma once

#include "includes/model_part.h"
#include "processes/process.h"

#include <unordered_map>

namespace Kratos
{

using NodeIdsToBoundariesHashMap     = std::unordered_multimap<std::vector<std::size_t>,
                                                               std::vector<GeometricalObject::Pointer>,
                                                               KeyHasherRange<std::vector<std::size_t>>,
                                                               KeyComparorRange<std::vector<std::size_t>>>;
using SortedToUnsortedNodeIdsHashMap = std::unordered_multimap<std::vector<std::size_t>,
                                                               std::vector<std::size_t>,
                                                               KeyHasherRange<std::vector<std::size_t>>,
                                                               KeyComparorRange<std::vector<std::size_t>>>;

class KRATOS_API(GEO_MECHANICS_APPLICATION) FindNeighbourElementsOfConditionsProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(FindNeighbourElementsOfConditionsProcess);

    explicit FindNeighbourElementsOfConditionsProcess(ModelPart& rModelPart);
    FindNeighbourElementsOfConditionsProcess& operator=(const FindNeighbourElementsOfConditionsProcess&) = delete;
    FindNeighbourElementsOfConditionsProcess(const FindNeighbourElementsOfConditionsProcess&) = delete;
    ~FindNeighbourElementsOfConditionsProcess() override = default;

    void                      Execute() override;
    void                      ExecuteInitialize() override;
    [[nodiscard]] std::string Info() const override;
    void                      PrintData(std::ostream& rOStream) const override;

private:
    ModelPart&                     mrModelPart;
    NodeIdsToBoundariesHashMap     mBoundaryNodeIdsToBoundaries;
    SortedToUnsortedNodeIdsHashMap mSortedToUnsortedBoundaryNodeIds;

    void InitializeBoundaryMaps();
    void FindNeighbouringElementsForAllBoundaryTypes();

    void SetElementAsNeighbourOfAllGeometryWithIdenticalNodeIds(const std::vector<std::size_t>& rNodeIds,
                                                                Element* pElement);

    void FindNeighboursBasedOnBoundaryType(auto GenerateBoundaries);

    void AddNeighbouringElementsBasedOnOverlappingBoundaryGeometries(
        Element& rElement, const Geometry<Node>::GeometriesArrayType& rElementBoundaryGeometries);

    void                      SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(Element& rElement,
                                                                                 const std::vector<std::size_t>& element_boundary_node_ids,
                                                                                 const GeometryData::KratosGeometryOrderType& r_order_type,
                                                                                 std::size_t LocalSpaceDimension);
    [[nodiscard]] static bool AreRotatedEquivalents(const std::vector<std::size_t>& rFirst,
                                                    const std::vector<std::size_t>& rSecond,
                                                    const GeometryData::KratosGeometryOrderType& rOrderType,
                                                    std::size_t LocalSpaceDimension);
    [[nodiscard]] static bool AreLinearRotatedEquivalents(std::vector<std::size_t> elements_boundary_node_ids,
                                                          const std::vector<std::size_t>& condition_node_ids);
    [[nodiscard]] static bool AreQuadraticRotatedEquivalents(std::vector<std::size_t> elements_boundary_node_ids,
                                                             const std::vector<std::size_t>& condition_node_ids,
                                                             std::size_t LocalSpaceDimension);

    [[nodiscard]] bool AllBoundariesHaveAtLeastOneNeighbour() const;
    [[nodiscard]] static std::vector<std::size_t> GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry);
    [[noreturn]] void ReportBoundariesWithoutNeighboursAndThrow() const;
};

std::ostream& operator<<(std::ostream& rOStream, const FindNeighbourElementsOfConditionsProcess& rThis);

} // namespace Kratos