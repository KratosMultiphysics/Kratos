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

#include "find_neighbour_elements_of_conditions_process.h"
#include "geometries/geometry.h"
#include <algorithm>
#include <iterator>

namespace Kratos
{

FindNeighbourElementsOfConditionsProcess::FindNeighbourElementsOfConditionsProcess(ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

void FindNeighbourElementsOfConditionsProcess::Execute()
{
    if (mrModelPart.Conditions().empty()) return;

    InitializeConditionMaps();
    FindNeighbouringElementsForAllBoundaryTypes();

    if (!AllConditionsHaveAtLeastOneNeighbour()) {
        ReportConditionsWithoutNeighboursAndThrow();
    }
}

void FindNeighbourElementsOfConditionsProcess::InitializeConditionMaps()
{
    mConditionNodeIdsToConditions.clear();
    std::ranges::transform(
        mrModelPart.Conditions(),
        std::inserter(mConditionNodeIdsToConditions, mConditionNodeIdsToConditions.end()), [](auto& rCondition) {
        return NodeIdsToConditionsHashMap::value_type(
            GetNodeIdsFromGeometry(rCondition.GetGeometry()), {&rCondition});
    });

    mSortedToUnsortedConditionNodeIds.clear();
    std::ranges::transform(
        mConditionNodeIdsToConditions,
        std::inserter(mSortedToUnsortedConditionNodeIds, mSortedToUnsortedConditionNodeIds.end()),
        [](const auto& rPair) {
        auto sorted_ids = rPair.first;
        std::ranges::sort(sorted_ids);
        return std::make_pair(sorted_ids, rPair.first);
    });
}

void FindNeighbourElementsOfConditionsProcess::FindNeighbouringElementsForAllBoundaryTypes()
{
    auto generate_generic_boundaries = [](const auto& rGeometry) {
        return rGeometry.GenerateBoundariesEntities();
    };
    auto generate_points   = [](const auto& rGeometry) { return rGeometry.GeneratePoints(); };
    auto generate_edges_3d = [](const auto& rGeometry) {
        return rGeometry.LocalSpaceDimension() == 3 ? rGeometry.GenerateEdges()
                                                    : PointerVector<Geometry<Node>>();
    };
    auto generate_edges_1d = [](const auto& rGeometry) {
        return rGeometry.LocalSpaceDimension() == 1 ? rGeometry.GenerateEdges()
                                                    : PointerVector<Geometry<Node>>();
    };

    // Note the order in the generators: the 1D elements are only added
    // as neighbours when the condition is not neighbouring 2D or 3D elements
    const std::vector<std::function<PointerVector<Geometry<Node>>(const Geometry<Node>&)>> boundary_generators = {
        generate_generic_boundaries, generate_points, generate_edges_3d, generate_edges_1d};

    for (const auto& r_boundary_generator : boundary_generators) {
        FindConditionNeighboursBasedOnBoundaryType(r_boundary_generator);
        if (AllConditionsHaveAtLeastOneNeighbour()) return;
    }
}

void FindNeighbourElementsOfConditionsProcess::FindConditionNeighboursBasedOnBoundaryType(auto GenerateBoundaries)
{
    for (auto& r_element : mrModelPart.Elements()) {
        const auto& r_element_geometry  = r_element.GetGeometry();
        const auto  boundary_geometries = GenerateBoundaries(r_element_geometry);
        AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(r_element, boundary_geometries);
    }
}

void FindNeighbourElementsOfConditionsProcess::AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
    Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries)
{
    for (const auto& r_boundary_geometry : rBoundaryGeometries) {
        const auto element_boundary_node_ids = GetNodeIdsFromGeometry(r_boundary_geometry);

        if (mConditionNodeIdsToConditions.contains(element_boundary_node_ids)) {
            SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(element_boundary_node_ids, &rElement);
        } else if (r_boundary_geometry.LocalSpaceDimension() == 2) {
            // No condition is directly found for this boundary, but it might be a rotated equivalent
            SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
                rElement, element_boundary_node_ids, r_boundary_geometry.GetGeometryOrderType());
        }
    }
}

void FindNeighbourElementsOfConditionsProcess::SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(
    const std::vector<std::size_t>& rConditionNodeIds, Element* pElement)
{
    const auto [start, end] = mConditionNodeIdsToConditions.equal_range(rConditionNodeIds);
    for (auto it = start; it != end; ++it) {
        const auto& r_conditions  = it->second;
        auto vector_of_neighbours = GlobalPointersVector<Element>{Element::WeakPointer{pElement}};

        for (auto& p_condition : r_conditions) {
            p_condition->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
        }
    }
}

void FindNeighbourElementsOfConditionsProcess::SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
    Element&                                     rElement,
    const std::vector<std::size_t>&              element_boundary_node_ids,
    const GeometryData::KratosGeometryOrderType& r_order_type)
{
    auto sorted_boundary_node_ids = element_boundary_node_ids;
    std::ranges::sort(sorted_boundary_node_ids);

    if (mSortedToUnsortedConditionNodeIds.contains(sorted_boundary_node_ids)) {
        const auto unsorted_condition_node_ids =
            mSortedToUnsortedConditionNodeIds.find(sorted_boundary_node_ids)->second;
        if (AreRotatedEquivalents(element_boundary_node_ids, unsorted_condition_node_ids, r_order_type)) {
            SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(unsorted_condition_node_ids, &rElement);
        }
    }
}

bool FindNeighbourElementsOfConditionsProcess::AreRotatedEquivalents(const std::vector<std::size_t>& rFirst,
                                                                     const std::vector<std::size_t>& rSecond,
                                                                     const GeometryData::KratosGeometryOrderType& rOrderType)
{
    switch (rOrderType) {
        using enum GeometryData::KratosGeometryOrderType;
    case Kratos_Linear_Order:
        return AreLinearRotatedEquivalents(rFirst, rSecond);
    case Kratos_Quadratic_Order:
        return AreQuadraticRotatedEquivalents(rFirst, rSecond);
    default:
        return false;
    }
}

bool FindNeighbourElementsOfConditionsProcess::AreLinearRotatedEquivalents(std::vector<std::size_t> First,
                                                                           const std::vector<std::size_t>& rSecond)
{
    const auto amount_of_needed_rotations = std::ranges::find(First, rSecond[0]) - First.begin();
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, First.end());
    return First == rSecond;
}

bool FindNeighbourElementsOfConditionsProcess::AreQuadraticRotatedEquivalents(std::vector<std::size_t> First,
                                                                              const std::vector<std::size_t>& rSecond)
{
    const auto amount_of_needed_rotations = std::ranges::find(First, rSecond[0]) - First.begin();
    auto       first_mid_side_node_id     = First.begin() + First.size() / 2;
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, first_mid_side_node_id);

    std::rotate(first_mid_side_node_id, first_mid_side_node_id + amount_of_needed_rotations, First.end());

    return First == rSecond;
}

bool FindNeighbourElementsOfConditionsProcess::AllConditionsHaveAtLeastOneNeighbour() const
{
    return std::ranges::all_of(mrModelPart.Conditions(), [](auto& rCondition) {
        return rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() > 0;
    });
}

void FindNeighbourElementsOfConditionsProcess::ReportConditionsWithoutNeighboursAndThrow() const
{
    std::vector<Condition> conditions_without_neighbours;
    std::ranges::copy_if(
        mrModelPart.Conditions(), std::back_inserter(conditions_without_neighbours),
        [](const auto& rCondition) { return rCondition.GetValue(NEIGHBOUR_ELEMENTS).size() == 0; });

    std::vector<std::size_t> ids_of_conditions_without_neighbours;
    std::ranges::transform(conditions_without_neighbours,
                           std::back_inserter(ids_of_conditions_without_neighbours),
                           [](const auto& rCondition) { return rCondition.Id(); });

    KRATOS_ERROR << "The condition(s) with the following ID(s) is/are found without any "
                    "corresponding element: "
                 << ids_of_conditions_without_neighbours << std::endl;
}

std::vector<std::size_t> FindNeighbourElementsOfConditionsProcess::GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry)
{
    std::vector<std::size_t> result;
    result.reserve(rGeometry.size());
    std::ranges::transform(rGeometry, std::back_inserter(result),
                           [](const auto& rNode) { return rNode.Id(); });
    return result;
}

std::ostream& operator<<(std::ostream& rOStream, const FindNeighbourElementsOfConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

std::string FindNeighbourElementsOfConditionsProcess::Info() const
{
    return "FindNeighbourElementsOfConditionsProcess";
}

void FindNeighbourElementsOfConditionsProcess::PrintData(std::ostream& rOStream) const
{
    this->PrintInfo(rOStream);
}

} // namespace Kratos