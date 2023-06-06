// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes



// Project includes
#include "includes/model_part.h"
#include "custom_processes/element_deactivation_process.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

ElementDeactivationProcess::ElementDeactivationProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ) : mrThisModelPart(rThisModelPart),
    mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mVariableName = mThisParameters["variable_name"].GetString();
    mThreshold = mThisParameters["variable_maximum_threshold"].GetDouble();
    mAverageOverIP = mThisParameters["average_calculation_over_ip"].GetBool();
}

/***********************************************************************************/
/***********************************************************************************/

void ElementDeactivationProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    if (KratosComponents<Variable<double>>::Has(mVariableName)) {
        // double type variable
        const auto &r_variable = KratosComponents<Variable<double>>::Get(mVariableName);

        block_for_each(mrThisModelPart.Elements(), [&](Element& rElement) {
            std::vector<double> element_data;
            rElement.CalculateOnIntegrationPoints(r_variable, element_data, mrThisModelPart.GetProcessInfo());

            if (mAverageOverIP) {
                double average_value = element_data[0];
                for (IndexType i = 1; i < element_data.size(); ++i)
                    average_value += element_data[i];
                
                if (average_value >= mThreshold)
                    rElement.Set(ACTIVE, false);
            } else {
                IndexType failed_ip_counter = 0;
                for (IndexType i = 0; i < element_data.size(); ++i) {
                    if (element_data[i] >= mThreshold)
                        failed_ip_counter++;
                }
                if (failed_ip_counter == element_data.size())
                    rElement.Set(ACTIVE, false);
            }
        });
    } else if (KratosComponents<Variable<Vector>>::Has(mVariableName)) {

    }






    KRATOS_CATCH("")

}

/***********************************************************************************/
/***********************************************************************************/

const Parameters ElementDeactivationProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"             : "please_specify_model_part_name",
        "variable_name"               : "DAMAGE",
        "variable_maximum_threshold"  : 0.9999,
        "average_calculation_over_ip" : true
    })");

    return default_parameters;
}

} // namespace Kratos.