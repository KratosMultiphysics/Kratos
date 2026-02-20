//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolo Crescenzio

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Application includes
#include "mpm_application_variables.h"
#include "reaction_utilities.h"

namespace Kratos
{
    array_1d<double,3> ReactionUtilities::CalculateGridConformingReaction(ModelPart& rModelPart) {
        return VariableUtils().SumHistoricalVariable<array_1d<double,3>>(REACTION, rModelPart, 0);
    }

    array_1d<double,3> ReactionUtilities::CalculateNonConformingReaction(ModelPart& rModelPart) {
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        std::vector<array_1d<double,3>> mpc_reaction{ ZeroVector(3) };
        return block_for_each<SumReduction<array_1d<double,3>>>(rModelPart.Conditions(), mpc_reaction,
            [&r_current_process_info](auto& condition, auto& rReaction) {
                condition.CalculateOnIntegrationPoints(MPC_CONTACT_FORCE, rReaction, r_current_process_info);
                return rReaction[0];
            });
    }
}
