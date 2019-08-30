//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//

// Project includes
#include "check_wake_condition_process.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos {

template <int Dim>
void CheckWakeConditionProcess<Dim>::Execute()
{
    unsigned int number_of_unfulfilled_wake_conditions = 0;
    for (auto& r_element : mrWakeModelPart.Elements()){
        const bool wake_condition_is_fulfilled =
            PotentialFlowUtilities::CheckIfWakeConditionIsFulfilled<Dim, Dim+1>(
                r_element, mTolerance, mEchoLevel);
        if (!wake_condition_is_fulfilled){
            number_of_unfulfilled_wake_conditions += 1;
        }
    }
    KRATOS_WARNING_IF("CheckIfWakeConditionIsFulfilled", number_of_unfulfilled_wake_conditions > 0)
            << "THE WAKE CONDITION IS NOT FULFILLED IN " << number_of_unfulfilled_wake_conditions << " ELEMENTS" << std::endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class CheckWakeConditionProcess<2>;
template class CheckWakeConditionProcess<3>;
} // namespace Kratos.
