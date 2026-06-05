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
//                   Richard Faasse
//

#include "find_neighbours_of_interfaces_process.h"
#include "containers/model.h"
#include "custom_utilities/neighbouring_element_finder.hpp"
#include "custom_utilities/process_utilities.h"
#include "includes/kratos_parameters.h"

using namespace std::string_literals;

namespace Kratos
{
FindNeighboursOfInterfacesProcess::FindNeighboursOfInterfacesProcess(Model& rModel, const Parameters& rProcessSettings)
    : mrInterfaceModelParts{ProcessUtilities::GetModelPartsFromSettings(
          rModel, rProcessSettings, FindNeighboursOfInterfacesProcess::Info())},
      mrModelPartWithCandidateNeighbours(
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

    for (const auto& r_interface_model_part : mrInterfaceModelParts) {
        element_finder.FindEntityNeighbours(r_interface_model_part.get().Elements(),
                                            mrModelPartWithCandidateNeighbours.Elements(),
                                            boundary_generator_by_local_dim);
    }
}

void FindNeighboursOfInterfacesProcess::RemoveNeighboursWithoutHigherLocalDimension() const
{
    for (const auto& r_interface_model_part : mrInterfaceModelParts) {
        for (auto& r_interface_element : r_interface_model_part.get().Elements()) {
            auto&      r_neighbour_elements = r_interface_element.GetValue(NEIGHBOUR_ELEMENTS);
            const auto interface_element_local_dimension =
                r_interface_element.GetGeometry().LocalSpaceDimension();
            auto is_neighbour_without_higher_local_dimension =
                [interface_element_local_dimension](const Element& rNeighbourElement) {
                return rNeighbourElement.GetGeometry().LocalSpaceDimension() <= interface_element_local_dimension;
            };
            // WORKAROUND: std::erase-remove on GlobalPointersVector corrupts adjacent memory.
            // Instead of erasing in-place, build a new container with elements to keep.
            GlobalPointersVector<Element> temp_neighbours;
            for (auto it = r_neighbour_elements.ptr_begin(); it != r_neighbour_elements.ptr_end(); ++it) {
                if (!is_neighbour_without_higher_local_dimension(**it)) {
                    temp_neighbours.push_back(*it);
                }
            }
            r_neighbour_elements = temp_neighbours;
        }
    }
}

void FindNeighboursOfInterfacesProcess::ExecuteFinalize()
{
    for (const auto& r_interface_model_part : mrInterfaceModelParts) {
        for (auto& r_interface_element : r_interface_model_part.get().Elements()) {
            r_interface_element.GetValue(NEIGHBOUR_ELEMENTS).clear();
        }
    }
}

std::string FindNeighboursOfInterfacesProcess::Info() const
{
    return "FindNeighboursOfInterfacesProcess"s;
}

} // namespace Kratos
