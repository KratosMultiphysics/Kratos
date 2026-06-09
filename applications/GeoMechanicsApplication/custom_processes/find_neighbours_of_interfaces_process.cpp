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
namespace
{
GlobalPointersVector<Element> KeepNeighboursWithHigherLocalDimension(const GlobalPointersVector<Element>& rNeighbourElements,
                                                                     const std::size_t InterfaceElementLocalDimension)
{
    GlobalPointersVector<Element> kept_neighbours;
    for (auto it = rNeighbourElements.ptr_begin(); it != rNeighbourElements.ptr_end(); ++it) {
        if ((**it).GetGeometry().LocalSpaceDimension() > InterfaceElementLocalDimension) {
            kept_neighbours.push_back(*it);
        }
    }

    return kept_neighbours;
}
} // namespace

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
            // WORKAROUND: std::erase-remove on GlobalPointersVector corrupts adjacent memory.
            // Instead of erasing in-place, construct a filtered container and assign.
            r_neighbour_elements = KeepNeighboursWithHigherLocalDimension(
                r_neighbour_elements, interface_element_local_dimension);
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
