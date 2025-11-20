// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "find_neighbours_of_interfaces_process.h"
#include "containers/model.h"
#include "custom_utilities/neighbouring_element_finder.hpp"
#include "custom_utilities/process_utilities.h"
#include "includes/kratos_parameters.h"

#include <string>

using namespace std::string_literals;

namespace Kratos
{
FindNeighboursOfInterfacesProcess::FindNeighboursOfInterfacesProcess(Model& rModel, const Parameters& rProcessSettings)
    : mrModelParts{ProcessUtilities::GetModelPartsFromSettings(
          rModel, rProcessSettings, FindNeighboursOfInterfacesProcess::Info())},
      mrMainModelPart(
          rModel.GetModelPart(rProcessSettings["model_part_name_for_neighbouring_elements"].GetString()))
{
}

void FindNeighboursOfInterfacesProcess::ExecuteInitialize()
{
    FindAllNeighboursOfElements();
    FilterOutNeighboursWhichDoNotHaveHigherLocalDimension();
}

void FindNeighboursOfInterfacesProcess::FindAllNeighboursOfElements()
{
    constexpr auto                                         enable_reverse_search = true;
    NeighbouringElementFinder                              element_finder(enable_reverse_search);
    NeighbouringElementFinder::BoundaryGeneratorByLocalDim boundary_generator_by_local_dim;
    boundary_generator_by_local_dim[std::size_t{1}] = std::make_unique<EdgesGenerator>();

    for (const auto& r_model_part : mrModelParts) {
        element_finder.FindEntityNeighbours(r_model_part.get().Elements(), mrMainModelPart.Elements(),
                                            boundary_generator_by_local_dim);
    }
}

void FindNeighboursOfInterfacesProcess::FilterOutNeighboursWhichDoNotHaveHigherLocalDimension() const
{
    for (const auto& r_model_part : mrModelParts) {
        for (auto& r_element : r_model_part.get().Elements()) {
            auto& r_neighbour_elements = r_element.GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.erase(
                std::remove_if(r_neighbour_elements.begin(), r_neighbour_elements.end(),
                               [&r_element](const Element& rNeighbourElement) {
                return rNeighbourElement.GetGeometry().LocalSpaceDimension() <=
                       r_element.GetGeometry().LocalSpaceDimension();
            }),
                r_neighbour_elements.end());
        }
    }
}

void FindNeighboursOfInterfacesProcess::ExecuteFinalize()
{
    for (const auto& r_model_part : mrModelParts) {
        for (auto& r_element : r_model_part.get().Elements()) {
            r_element.GetValue(NEIGHBOUR_ELEMENTS).clear();
        }
    }
}

std::string FindNeighboursOfInterfacesProcess::Info() const
{
    return "FindNeighboursOfInterfacesProcess"s;
}

} // namespace Kratos
