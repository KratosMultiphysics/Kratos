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
//

#include "custom_processes/deactivate_conditions_on_inactive_elements_process.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
using namespace std::string_literals;

DeactivateConditionsOnInactiveElements::DeactivateConditionsOnInactiveElements(ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart)
{
}

void DeactivateConditionsOnInactiveElements::Execute()
{
    auto is_active = [](const auto& rpNeighbourElement) {
        return rpNeighbourElement->IsDefined(ACTIVE) ? rpNeighbourElement->Is(ACTIVE) : true;
    };

    block_for_each(mrModelPart.Conditions(), [&is_active](Condition& rCondition) {
        const auto& vector_of_neighbours = rCondition.GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_ERROR_IF(vector_of_neighbours.size() == 0)
            << "Condition without any corresponding element, ID " << rCondition.Id() << "\n"
            << "Call a process to find neighbour elements before calling this function." << std::endl;

        rCondition.Set(ACTIVE, std::any_of(vector_of_neighbours.ptr_begin(),
                                           vector_of_neighbours.ptr_end(), is_active));
    });
}

std::string DeactivateConditionsOnInactiveElements::Info() const
{
    return "DeactivateConditionsOnInactiveElements"s;
}

void DeactivateConditionsOnInactiveElements::PrintData(std::ostream& rOStream) const
{
    this->PrintInfo(rOStream);
}

} // namespace Kratos