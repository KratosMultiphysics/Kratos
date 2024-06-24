// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors: Aron Noordam,
//
#include "calculate_incremental_displacement_process.h"
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

    CalculateIncrementalDisplacementProcess::CalculateIncrementalDisplacementProcess(ModelPart& rModelPart, const Parameters& rProcessSettings)
    : Process(Flags()), mrModelPart{rModelPart}
{
}

void CalculateIncrementalDisplacementProcess::Execute()
{
    KRATOS_TRY

		block_for_each(mrModelPart.Nodes(), [&](Node& rNode) {
		rNode.FastGetSolutionStepValue(INCREMENTAL_DISPLACEMENT) = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0) - rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
			});
  
    KRATOS_CATCH("")
}

CalculateIncrementalDisplacementProcess::~CalculateIncrementalDisplacementProcess() = default;

} // namespace Kratos
