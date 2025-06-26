// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//                   Jonathan Nuttall
//

#include "fix_water_pressures_above_phreatic_line.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

namespace Kratos
{

FixWaterPressuresAbovePhreaticLineProcess::FixWaterPressuresAbovePhreaticLineProcess(ModelPart& rMainModelPart,
                                                                                     const Parameters& rSettings)
    : mrModelPart(rMainModelPart)
{
    const auto x_coordinates = rSettings["x_coordinates"].GetVector();
    const auto y_coordinates = rSettings["y_coordinates"].GetVector();

    for (int i = 0; i < x_coordinates.size(); ++i) {
        mPhreaticLineTable.insert(x_coordinates[i], y_coordinates[i]);
    }
}

void FixWaterPressuresAbovePhreaticLineProcess::ExecuteInitializeSolutionStep()
{
    block_for_each(mrModelPart.Nodes(), [this](Node& rNode) {
        const auto current_coordinates =
            rNode.GetInitialPosition() + rNode.FastGetSolutionStepValue(TOTAL_DISPLACEMENT);
        if (current_coordinates[1] > mPhreaticLineTable(current_coordinates[0])) {
            rNode.FastGetSolutionStepValue(WATER_PRESSURE) =
                0.0; // Could be changed to small positive value instead of zero
            rNode.Fix(WATER_PRESSURE);
        } else {
            rNode.Free(WATER_PRESSURE); // We could think about freeing in finalizesolutionstep
        }
    });
}

} // namespace Kratos