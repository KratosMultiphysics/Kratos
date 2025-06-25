//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Richard Faasse
//

#include "fix_water_pressures_above_phreatic_line.h"
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
        if (rNode.Y() > mPhreaticLineTable(rNode.X())) {
            rNode.FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
            rNode.Fix(WATER_PRESSURE);
        }
        else {
            rNode.Free(WATER_PRESSURE);
        }
    });
}

} // namespace Kratos