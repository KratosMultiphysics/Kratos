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

namespace Kratos
{

using NodeIdToConditionsHashMap      = std::unordered_multimap<std::vector<std::size_t>,
                                                               std::vector<Condition::Pointer>,
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

    explicit FindNeighbourElementsOfConditionsProcess(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    FindNeighbourElementsOfConditionsProcess& operator=(const FindNeighbourElementsOfConditionsProcess&) = delete;
    FindNeighbourElementsOfConditionsProcess(const FindNeighbourElementsOfConditionsProcess&) = delete;
    ~FindNeighbourElementsOfConditionsProcess() override = default;

    void Execute() override;

    [[nodiscard]] std::string Info() const override
    {
        return "FindNeighbourElementsOfConditionsProcess";
    }

    void PrintData(std::ostream& rOStream) const override { this->PrintInfo(rOStream); }

private:
    ModelPart&                     mrModelPart;
    NodeIdToConditionsHashMap      mConditionNodeIdsToCondition;
    SortedToUnsortedNodeIdsHashMap mSortedToUnsortedConditionNodeIds;

    void InitializeConditionMaps();
    void FindNeighbouringElementsForAllBoundaryTypes();

    void SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(const std::vector<std::size_t>& rConditionNodeIds,
                                                                  Element* pElement);

    void FindConditionNeighboursBasedOnBoundaryType(auto generate_boundaries);

    void AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
        Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries);

    void                      SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(Element& rElement,
                                                                                 const std::vector<std::size_t>& element_boundary_node_ids,
                                                                                 const GeometryData::KratosGeometryOrderType& r_order_type);
    [[nodiscard]] static bool AreRotatedEquivalents(const std::vector<std::size_t>& rFirst,
                                                    const std::vector<std::size_t>& rSecond,
                                                    const GeometryData::KratosGeometryOrderType& rOrderType);
    [[nodiscard]] static bool AreLinearRotatedEquivalents(std::vector<std::size_t> elements_boundary_node_ids,
                                                          const std::vector<std::size_t>& condition_node_ids);
    [[nodiscard]] static bool AreQuadraticRotatedEquivalents(std::vector<std::size_t> elements_boundary_node_ids,
                                                             const std::vector<std::size_t>& condition_node_ids);

    [[nodiscard]] bool AllConditionsHaveAtLeastOneNeighbour() const;
    [[nodiscard]] static std::vector<std::size_t> GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry);
    void ReportConditionsWithoutNeighbours() const;
};

inline std::ostream& operator<<(std::ostream& rOStream, const FindNeighbourElementsOfConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos