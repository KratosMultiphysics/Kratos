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
#include "custom_utilities/neighbouring_element_finder.hpp"
#include "geometries/geometry.h"
#include "includes/model_part.h"

#include <algorithm>
#include <iterator>

namespace Kratos
{
using namespace std::string_literals;

FindNeighbourElementsOfConditionsProcess::FindNeighbourElementsOfConditionsProcess(ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

void FindNeighbourElementsOfConditionsProcess::Execute()
{
    if (mrModelPart.Conditions().empty()) return;

    FindNeighbouringElementsForAllBoundaryTypes();

    if (!AllConditionsHaveAtLeastOneNeighbour()) {
        ReportConditionsWithoutNeighboursAndThrow();
    }
}

void FindNeighbourElementsOfConditionsProcess::FindNeighbouringElementsForAllBoundaryTypes()
{
    NeighbouringElementFinder::BoundaryGeneratorByLocalDim boundary_generator_map;
    boundary_generator_map[std::size_t{0}] = std::make_unique<PointsGenerator>();
    boundary_generator_map[std::size_t{1}] = std::make_unique<EdgesGenerator>();
    boundary_generator_map[std::size_t{2}] = std::make_unique<FacesGenerator>();

    NeighbouringElementFinder finder;
    finder.FindEntityNeighbours(mrModelPart.Conditions(), mrModelPart.Elements(), boundary_generator_map);
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
    return "FindNeighbourElementsOfConditionsProcess"s;
}

void FindNeighbourElementsOfConditionsProcess::PrintData(std::ostream& rOStream) const
{
    this->PrintInfo(rOStream);
}

} // namespace Kratos