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

    using IndexType    = std::size_t;
    using NodeType     = Node;
    using GeometryType = Geometry<NodeType>;

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
    ModelPart& mrModelPart;

    [[nodiscard]] bool AllConditionsAreVisited() const;

    static void SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(NodeIdToConditionsHashMap&      rFacesMap,
                                                    const std::vector<std::size_t>& rConditionNodeIds,
                                                    Element*                        pElement);

    void AddNeighboringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
        NodeIdToConditionsHashMap&                 FacesMap,
        const SortedToUnsortedNodeIdsHashMap&      FacesMapSorted,
        Element&                                   rElement,
        const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries) const;

    [[nodiscard]] bool FindPermutations(std::vector<std::size_t>        elements_boundary_node_ids,
                                        const std::vector<std::size_t>& condition_node_ids) const;
    [[nodiscard]] bool FindPermutationsQuadratic(std::vector<std::size_t> elements_boundary_node_ids,
                                                 const std::vector<std::size_t>& condition_node_ids) const;
    void CheckBoundaryTypeForAllElements(auto                       generate_boundaries,
                                         NodeIdToConditionsHashMap& condition_node_ids_to_condition,
                                         SortedToUnsortedNodeIdsHashMap& sorted_condition_node_ids_to_condition);

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