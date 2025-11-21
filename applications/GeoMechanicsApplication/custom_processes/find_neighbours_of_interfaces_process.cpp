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
    RemoveNeighboursWithoutHigherLocalDimension();
}

void FindNeighboursOfInterfacesProcess::FindAllNeighboursOfElements()
{
    constexpr auto                                         enable_reverse_search = true;
    NeighbouringElementFinder                              element_finder(enable_reverse_search);
    NeighbouringElementFinder::BoundaryGeneratorByLocalDim boundary_generator_by_local_dim;
    boundary_generator_by_local_dim[std::size_t{1}] = std::make_unique<EdgesGenerator>();
    boundary_generator_by_local_dim[std::size_t{2}] = std::make_unique<FacesGenerator>();

    for (const auto& r_model_part : mrModelParts) {
        element_finder.FindEntityNeighbours(r_model_part.get().Elements(), mrMainModelPart.Elements(),
                                            boundary_generator_by_local_dim);
    }
}

void FindNeighboursOfInterfacesProcess::RemoveNeighboursWithoutHigherLocalDimension() const
{
    for (const auto& r_model_part : mrModelParts) {
        for (auto& r_element : r_model_part.get().Elements()) {
            auto&      r_neighbour_elements    = r_element.GetValue(NEIGHBOUR_ELEMENTS);
            const auto element_local_dimension = r_element.GetGeometry().LocalSpaceDimension();
            auto       is_neighbour_without_higher_local_dimension =
                [element_local_dimension](const Element& rNeighbourElement) {
                return rNeighbourElement.GetGeometry().LocalSpaceDimension() <= element_local_dimension;
            };
            r_neighbour_elements.erase(
                std::remove_if(r_neighbour_elements.begin(), r_neighbour_elements.end(),
                               is_neighbour_without_higher_local_dimension),
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
