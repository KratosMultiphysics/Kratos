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

    mNeighbouringEntityFinder.InitializeBoundaryMaps(mrModelPart.Conditions());
    FindNeighbouringElementsForAllBoundaryTypes();

    if (!AllConditionsHaveAtLeastOneNeighbour()) {
        ReportConditionsWithoutNeighboursAndThrow();
    }
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
        mNeighbouringEntityFinder.FindConditionNeighboursBasedOnBoundaryType(
            r_boundary_generator, mrModelPart.Elements());
        if (AllConditionsHaveAtLeastOneNeighbour()) return;
    }
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