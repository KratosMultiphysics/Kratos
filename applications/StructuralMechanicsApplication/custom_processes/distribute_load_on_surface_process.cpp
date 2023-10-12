// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Armin Geiser
//

// System includes

// External includes

// Project includes
#include "distribute_load_on_surface_process.h"
#include "utilities/interval_utility.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

DistributeLoadOnSurfaceProcess::DistributeLoadOnSurfaceProcess(
    ModelPart& rModelPart,
    Parameters Settings)
    : mrModelPart(rModelPart),
      mParameters(Settings)
{
    Parameters default_parameters(R"(
        {
            "help"            : "This process distributes a load on surface load conditions belonging to a modelpart. The load is distributed according to the surface area.",
            "model_part_name" : "please_specify_model_part_name",
            "interval"        : [0.0, 1e30],
            "load"            : [1.0, 0.0, 0.0]
        }  )"
    );

    IntervalUtility interval_utility(mParameters);

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    KRATOS_ERROR_IF(mParameters["load"].GetVector().size() != 3) << "'load' has to be a vector of doubles with size 3!" << std::endl;
}

void DistributeLoadOnSurfaceProcess::ExecuteInitializeSolutionStep()
{
    const double current_time = mrModelPart.GetProcessInfo().GetValue(TIME);

    IntervalUtility interval_utility(mParameters);
    if (interval_utility.IsInInterval(current_time)) {
        // Calculate the total area
        const auto& r_conditions_array = mrModelPart.Conditions();
        const double total_area = block_for_each<SumReduction<double>>(r_conditions_array, [&](Condition& rCond) {
            return rCond.GetGeometry().Area();
        });

        // Compute the total area
        const double global_total_area = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(total_area);

        // Getting force by area
        const Vector force_by_area = mParameters["load"].GetVector() / global_total_area;

        // Assign on conditions
        block_for_each(r_conditions_array, [&](Condition& rCond) {
            rCond.SetValue(SURFACE_LOAD, force_by_area);
        });
    }
}

}  // namespace Kratos.
