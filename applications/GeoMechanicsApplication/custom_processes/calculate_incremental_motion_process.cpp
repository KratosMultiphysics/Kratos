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
#include "calculate_incremental_motion_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "includes/node.h"

namespace Kratos
{

CalculateIncrementalMotionProcess::CalculateIncrementalMotionProcess(ModelPart&        rModelPart,
                                                                     const Parameters& rParameters)
    : Process(Flags()), mrModelPart{rModelPart}
{
    mBaseVariableName = rParameters["variable_name"].GetString();
    if (mBaseVariableName == "DISPLACEMENT") {
        mResultsVariableName = "INCREMENTAL_DISPLACEMENT";
    } else if (mBaseVariableName == "ROTATION") {
        mResultsVariableName = "INCREMENTAL_ROTATION";
    } else {
        KRATOS_ERROR << "Invalid variable name: " << rParameters["variable_name"].GetString()
                     << ". Expected DISPLACEMENT or ROTATION." << std::endl;
    }
}

void CalculateIncrementalMotionProcess::Execute()
{
    if (KratosComponents<Variable<array_1d<double, 3>>>::Has(mResultsVariableName) &&
        KratosComponents<Variable<array_1d<double, 3>>>::Has(mBaseVariableName)) {
        const Variable<array_1d<double, 3>>& rResultVariable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mResultsVariableName);

        const Variable<array_1d<double, 3>>& rBaseVariable =
            KratosComponents<Variable<array_1d<double, 3>>>::Get(mBaseVariableName);

        block_for_each(mrModelPart.Nodes(), [rResultVariable, rBaseVariable](Node& rNode) {
            rNode.FastGetSolutionStepValue(rResultVariable) =
                rNode.FastGetSolutionStepValue(rBaseVariable, 0) -
                rNode.FastGetSolutionStepValue(rBaseVariable, 1);
        });
    } else {
        KRATOS_ERROR << "Variables " << mResultsVariableName << " and/or " << mBaseVariableName
                     << " not found in the model part." << std::endl;
    }
}

std::string CalculateIncrementalMotionProcess::Info() const
{
    return "CalculateIncrementalMotionProcess";
}

} // namespace Kratos
