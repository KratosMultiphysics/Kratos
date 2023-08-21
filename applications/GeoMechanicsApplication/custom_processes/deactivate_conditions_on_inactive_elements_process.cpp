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

#include "includes/kratos_flags.h"
#include "utilities/parallel_utilities.h"
#include "custom_processes/deactivate_conditions_on_inactive_elements_process.hpp"

namespace Kratos
{
void DeactivateConditionsOnInactiveElements::Execute()
{
    KRATOS_TRY

    block_for_each(mrModelPart.Conditions(), [&](Condition& rCondition) {
        const auto &VectorOfNeighbours = rCondition.GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_ERROR_IF(VectorOfNeighbours.size() == 0)
            << "Condition without any corresponding element, ID " << rCondition.Id() << "\n"
            << "Call a process to find neighbour elements before calling this function."
            << std::endl;

        bool IsElementActive = false;
        for (unsigned int i=0; i < VectorOfNeighbours.size(); ++i) {
            if (VectorOfNeighbours[i].IsDefined(ACTIVE)) {
                if (VectorOfNeighbours[i].Is(ACTIVE)) IsElementActive = true;
            } else {
                IsElementActive = true;
            }
        }
        rCondition.Set(ACTIVE, IsElementActive);
    });

    KRATOS_CATCH("")
}

}