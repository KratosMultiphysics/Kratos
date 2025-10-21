// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License: geo_mechanics_application/license.txt
//
//  Main authors: Aron Noordam,
//
#include "calculate_total_motion_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"

namespace Kratos
{

CalculateTotalMotionProcess::CalculateTotalMotionProcess(ModelPart& rModelPart, const Parameters& rParameters)
    : Process(Flags()), mrModelPart{rModelPart}
{
    if (rParameters["variable_name"].GetString() == "DISPLACEMENT") {
        mResultsVariableName     = "TOTAL_DISPLACEMENT";
        mIncrementalVariableName = "INCREMENTAL_DISPLACEMENT";
    } else if (rParameters["variable_name"].GetString() == "ROTATION") {
        mResultsVariableName     = "TOTAL_ROTATION";
        mIncrementalVariableName = "INCREMENTAL_ROTATION";

    } else {
        KRATOS_ERROR << "Invalid variable name: " << rParameters["variable_name"].GetString()
                     << ". Expected DISPLACEMENT or ROTATION." << std::endl;
    }
}

void CalculateTotalMotionProcess::Execute()
{
    if (KratosComponents<Variable<array_1d<double, 3>>>::Has(mResultsVariableName) &&
        KratosComponents<Variable<array_1d<double, 3>>>::Has(mIncrementalVariableName)) {
        const Variable<array_1d<double, 3>>& rResultVariable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mResultsVariableName);

        const Variable<array_1d<double, 3>>& rIncrementalVariable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mIncrementalVariableName);

        block_for_each(mrModelPart.Nodes(), [rResultVariable, rIncrementalVariable](Node& rNode) {
            rNode.FastGetSolutionStepValue(rResultVariable) +=
                rNode.FastGetSolutionStepValue(rIncrementalVariable);
        });
    } else {
        KRATOS_ERROR << "Variables " << mResultsVariableName << " and/or "
                     << mIncrementalVariableName << " not found in the model part." << std::endl;
    }
}

std::string CalculateTotalMotionProcess::Info() const { return "CalculateTotalMotionProcess"; }

} // namespace Kratos
