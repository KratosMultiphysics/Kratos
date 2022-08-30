//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Armin Geiser
//

// System includes

// External includes

// Project includes
#include "distribute_load_on_surface_process.h"
#include "utilities/interval_utility.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

DistributeLoadOnSurfaceProcess::DistributeLoadOnSurfaceProcess(ModelPart& rModelPart,
                                                                Parameters Settings)
                                                                : mrModelPart(rModelPart),
                                                                mParameters(Settings)
{
    Parameters default_parameters(R"(
        {
            "help"            : "This process distributes a load on surface load conditions belonging to a modelpart. The load is distributed according to the surface area.",
            "model_part_name" : "please_specify_model_part_name",
            "interval"        : [0.0, 1e30],
            "load"           : [1.0, 0.0, 0.0]
        }  )"
    );

    IntervalUtility interval_utility(mParameters);

    mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    KRATOS_ERROR_IF(mParameters["load"].GetVector().size() != 3) <<
        "'load' has to be a vector of doubles with size 3!" << std::endl;
}

void DistributeLoadOnSurfaceProcess::ExecuteInitializeSolutionStep()
{
    const double current_time = mrModelPart.GetProcessInfo().GetValue(TIME);

    IntervalUtility interval_utility(mParameters);
    if (interval_utility.IsInInterval(current_time)) {
        double total_area = 0.0;
        for (auto& r_cond : mrModelPart.Conditions()) {
            total_area += r_cond.GetGeometry().Area();
        }

        const double global_total_area = mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(total_area);

        Vector force_by_area = mParameters["load"].GetVector() / global_total_area;

        for (auto& r_cond : mrModelPart.Conditions()) {
            const double area = r_cond.GetGeometry().Area();
            r_cond.SetValue(SURFACE_LOAD, force_by_area * area);
        }
    }
}

}  // namespace Kratos.
