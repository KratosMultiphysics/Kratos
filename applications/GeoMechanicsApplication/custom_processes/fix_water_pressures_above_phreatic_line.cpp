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
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"

namespace Kratos
{

FixWaterPressuresAbovePhreaticLineProcess::FixWaterPressuresAbovePhreaticLineProcess(ModelPart& rMainModelPart,
                                                                                     const Parameters& rSettings)
    : mrModelPart(rMainModelPart), mMoveMesh(rSettings["move_mesh"].GetBool())
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
        const auto coordinates =
            mMoveMesh ? rNode.Coordinates()
                      : rNode.Coordinates() + rNode.FastGetSolutionStepValue(TOTAL_DISPLACEMENT);
        if (coordinates[1] > mPhreaticLineTable(coordinates[0])) {
            rNode.FastGetSolutionStepValue(WATER_PRESSURE) =
                0.0; // Could be changed to small positive value instead of zero
            rNode.Fix(WATER_PRESSURE);
        } else {
            rNode.Free(WATER_PRESSURE);
        }
    });
}

} // namespace Kratos