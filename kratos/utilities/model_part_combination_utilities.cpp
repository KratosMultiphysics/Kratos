//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "utilities/model_part_combination_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
void ModelPartCombinationUtilities::CombineModelParts(Parameters ThisParameters)
{
    // Ensuring parameters
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters ModelPartCombinationUtilities::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_parts_list"         : [],
        "combined_model_part_name" : "CombinedModelParts"
    })" );
    return default_parameters;
}

}  // namespace Kratos.
