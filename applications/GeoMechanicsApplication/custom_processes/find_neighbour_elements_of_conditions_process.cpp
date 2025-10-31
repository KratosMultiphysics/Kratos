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

void FindNeighbourElementsOfConditionsProcess::ExecuteInitialize()
{
    KRATOS_INFO("FindNeighbourElementsOfConditionsProcess::ExecuteInitialize()") << "For model part " << mrModelPart.Name() << std::endl;
    FindNeighbourElementsOfConditionsProcess::Execute();
}

void FindNeighbourElementsOfConditionsProcess::Execute()
{
    KRATOS_INFO("FindNeighbourElementsOfConditionsProcess::Execute") << "For model part " << mrModelPart.Name() << std::endl;
    if (mrModelPart.Conditions().empty()) return;

    InitializeBoundaryMaps();
    FindNeighbouringElementsForAllBoundaryTypes();

    if (!AllBoundariesHaveAtLeastOneNeighbour()) {
        ReportBoundariesWithoutNeighboursAndThrow();
    }
}

void FindNeighbourElementsOfConditionsProcess::InitializeBoundaryMaps()
{
    mBoundaryNodeIdsToBoundaries.clear();
    std::ranges::transform(
        mrModelPart.Conditions(),
        std::inserter(mBoundaryNodeIdsToBoundaries, mBoundaryNodeIdsToBoundaries.end()), [](auto& rBoundary) {
        // van de conditie wordt alleen de geometrie gebruikt
        return NodeIdsToBoundariesHashMap::value_type(
            GetNodeIdsFromGeometry(rBoundary.GetGeometry()), {&rBoundary});
    });

    mSortedToUnsortedBoundaryNodeIds.clear();
    std::ranges::transform(
        mBoundaryNodeIdsToBoundaries,
        std::inserter(mSortedToUnsortedBoundaryNodeIds, mSortedToUnsortedBoundaryNodeIds.end()),
        [](const auto& rPair) {
        auto sorted_ids = rPair.first;
        // nu wordt sorted_ids een oplopend lijstje, dat vergelijkt makkelijk. maar make_pair snap ik nog niet
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
    // as neighbours when the boundary is not neighbouring 2D or 3D elements
    const std::vector<std::function<PointerVector<Geometry<Node>>(const Geometry<Node>&)>> boundary_generators = {
        generate_generic_boundaries, generate_points, generate_edges_3d, generate_edges_1d};
    for (const auto& r_boundary_generator : boundary_generators) {
        FindNeighboursBasedOnBoundaryType(r_boundary_generator);
        if (AllBoundariesHaveAtLeastOneNeighbour()) return;
    }
}

void FindNeighbourElementsOfConditionsProcess::FindNeighboursBasedOnBoundaryType(auto GenerateBoundaries)
{
    for (auto& r_element : mrModelPart.Elements()) {
        const auto element_boundary_geometries = GenerateBoundaries(r_element.GetGeometry());
        AddNeighbouringElementsBasedOnOverlappingBoundaryGeometries(r_element, element_boundary_geometries);
    }
}

// dit lijkt dubbele input, via rElement denk ik bij de BoundaryGeometries van dat element te kunnen via GenerateBoundaries(r_element.GetGeometry())
// dan moet de GenerateBoundaries functie wel worden meegegeven
void FindNeighbourElementsOfConditionsProcess::AddNeighbouringElementsBasedOnOverlappingBoundaryGeometries(
    Element& rElement, const Geometry<Node>::GeometriesArrayType& rElementBoundaryGeometries)
{
    for (const auto& r_element_boundary_geometry: rElementBoundaryGeometries) {
        const auto element_boundary_node_ids = GetNodeIdsFromGeometry(r_element_boundary_geometry);

        if (mBoundaryNodeIdsToBoundaries.contains(element_boundary_node_ids)) {
            // the boundaries are accessed through the member mBoundaryNodeIdsToBoundaries
            SetElementAsNeighbourOfAllGeometryWithIdenticalNodeIds(element_boundary_node_ids, &rElement);
        } else // if (r_element_boundary_geometry.LocalSpaceDimension() == 2)
        {
            // No element boundary is directly found for this boundary, but it might be a rotated equivalent
            SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
                rElement, element_boundary_node_ids, r_element_boundary_geometry.GetGeometryOrderType(),
                r_element_boundary_geometry.LocalSpaceDimension());
        }
        // na toevoegen van else if (r_boundary_geometry.LocalSpaceDimension() == 1) { SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent()
        // vind je ook de omgekeerde edges van 2D/plane strain elementen
    }
}

// kan gegeneraliseerd worden rConditionNodeIds --> rGeometryNodeIds ( NB een geometry is al een geometry of nodes )
void FindNeighbourElementsOfConditionsProcess::SetElementAsNeighbourOfAllGeometryWithIdenticalNodeIds(
    const std::vector<std::size_t>& rNodeIds, Element* pElement)
{
    const auto [start, end] = mBoundaryNodeIdsToBoundaries.equal_range(rNodeIds);
    for (auto it = start; it != end; ++it) {
        const auto& r_boundaries = it->second;
        // dit snap ik niet
        auto vector_of_neighbours = GlobalPointersVector<Element>{Element::WeakPointer{pElement}};

        for (auto& p_boundary : r_boundaries) {
            // geeft SetValue een uitbreiding van het lijstje NEIGHBOUR_ELEMENTS of overschrijven we hier het geheel met een 1 lange vector?
            KRATOS_INFO("FindNeighbourElementsOfConditionsProcess::SetElementAsNeighbourOfAllGeometryWithIdenticalNodeIds") << vector_of_neighbours[0].Id() << std::endl;
            p_boundary->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
        }
    }
}

void FindNeighbourElementsOfConditionsProcess::SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
    Element&                                     rElement,
    const std::vector<std::size_t>&              element_boundary_node_ids,
    const GeometryData::KratosGeometryOrderType& r_order_type,
    std::size_t                                  LocalSpaceDimension)
{
    auto sorted_boundary_node_ids = element_boundary_node_ids;
    std::ranges::sort(sorted_boundary_node_ids);

    if (mSortedToUnsortedBoundaryNodeIds.contains(sorted_boundary_node_ids)) {
        const auto unsorted_boundary_node_ids =
            mSortedToUnsortedBoundaryNodeIds.find(sorted_boundary_node_ids)->second;
        if (AreRotatedEquivalents(element_boundary_node_ids, unsorted_boundary_node_ids,
                                  r_order_type, LocalSpaceDimension)) {
            SetElementAsNeighbourOfAllGeometryWithIdenticalNodeIds(unsorted_boundary_node_ids, &rElement);
        }
    }
}

bool FindNeighbourElementsOfConditionsProcess::AreRotatedEquivalents(const std::vector<std::size_t>& rFirst,
                                                                     const std::vector<std::size_t>& rSecond,
                                                                     const GeometryData::KratosGeometryOrderType& rOrderType,
                                                                     std::size_t LocalSpaceDimension)
{
    // finds if the geometry is a permutation. Works for plane and almost for line geometries. Currently only used for planes.
    switch (rOrderType) {
        using enum GeometryData::KratosGeometryOrderType;
    case Kratos_Linear_Order:
        return AreLinearRotatedEquivalents(rFirst, rSecond);
    case Kratos_Quadratic_Order:
        return AreQuadraticRotatedEquivalents(rFirst, rSecond, LocalSpaceDimension);
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

bool FindNeighbourElementsOfConditionsProcess::AreQuadraticRotatedEquivalents(
    std::vector<std::size_t> First, const std::vector<std::size_t>& rSecond, std::size_t LocalSpaceDimension)
{
    const auto amount_of_needed_rotations = std::ranges::find(First, rSecond[0]) - First.begin();
    auto first_mid_side_node_id = (LocalSpaceDimension == 1) ? First.begin()+2 : First.begin() + First.size() / 2;
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, first_mid_side_node_id);
    // Only rotate midside nodes of planes, the midside node for quadratic line remains in place
    if (LocalSpaceDimension == 2) {
        std::rotate(first_mid_side_node_id, first_mid_side_node_id + amount_of_needed_rotations,
                    First.end());
    }

    return First == rSecond;
}

bool FindNeighbourElementsOfConditionsProcess::AllBoundariesHaveAtLeastOneNeighbour() const
{
    return std::ranges::all_of(mrModelPart.Conditions(), [](auto& rBoundary) {
        return rBoundary.GetValue(NEIGHBOUR_ELEMENTS).size() > 0;
    });
}

void FindNeighbourElementsOfConditionsProcess::ReportBoundariesWithoutNeighboursAndThrow() const
{
    // i.p.v. mrModelPart.Conditions kan hier een lijst met geometry(Pointer?)s is als argument
    std::vector<GeometricalObject> boundaries_without_neighbours;
    std::ranges::copy_if(
        mrModelPart.Conditions(), std::back_inserter(boundaries_without_neighbours),
        [](const auto& rBoundary) { return rBoundary.GetValue(NEIGHBOUR_ELEMENTS).size() == 0; });

    std::vector<std::size_t> ids_of_boundaries_without_neighbours;
    std::ranges::transform(boundaries_without_neighbours,
                           std::back_inserter(ids_of_boundaries_without_neighbours),
                           [](const auto& rBoundary) { return rBoundary.Id(); });

    KRATOS_ERROR << "The condition(s) with the following ID(s) is/are found without any "
                    "corresponding element: "
                 << ids_of_boundaries_without_neighbours << std::endl;
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