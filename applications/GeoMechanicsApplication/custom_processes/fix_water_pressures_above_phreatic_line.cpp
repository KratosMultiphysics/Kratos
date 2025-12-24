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
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

#include <string>

namespace Kratos
{
using namespace std::string_literals;

FixWaterPressuresAbovePhreaticLineProcess::FixWaterPressuresAbovePhreaticLineProcess(ModelPart& rMainModelPart,
                                                                                     const Parameters& rSettings)
    : mrModelPart(rMainModelPart)
{
    const auto x_coordinates = rSettings["x_coordinates"].GetVector();
    const auto y_coordinates = rSettings["y_coordinates"].GetVector();

    KRATOS_ERROR_IF(x_coordinates.size() != y_coordinates.size())
        << "The lengths of the x_coordinates and y_coordinates of the phreatic line must be "
           "equal.\n";

    KRATOS_ERROR_IF(x_coordinates.empty())
        << "The x_coordinates and y_coordinates of the phreatic line must be "
           "non-empty vectors.\n";

    for (std::size_t i = 0; i < x_coordinates.size(); ++i) {
        mPhreaticLineTable.insert(x_coordinates[i], y_coordinates[i]);
    }
}

void FixWaterPressuresAbovePhreaticLineProcess::ExecuteInitializeSolutionStep()
{
    block_for_each(mrModelPart.Nodes(), [this](Node& rNode) {
        if (rNode.Y() >= mPhreaticLineTable.GetValue(rNode.X())) {
            rNode.Fix(WATER_PRESSURE);
            rNode.FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
        } else {
            rNode.Free(WATER_PRESSURE);
        }
    });
}

std::string FixWaterPressuresAbovePhreaticLineProcess::Info() const
{
    return "FixWaterPressuresAbovePhreaticLineProcess"s;
}

} // namespace Kratos