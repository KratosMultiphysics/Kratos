//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "response_utils.h"

namespace Kratos {

void ResponseUtils::CheckAndPrepareModelPartsForSensitivityComputation(
    const std::vector<ModelPart*>& rEvaluatedModelParts,
    const SensitivityModelPartVariablesListMap& rSensitivityModelPartVariableInfo,
    const Flags& rFlag,
    const std::vector<SensitivityFieldVariableTypes>& rUsedNodalSensitivityVariables)
{
    KRATOS_TRY

    // reset entity flags for sensitivity model parts
    for (auto& it : rSensitivityModelPartVariableInfo) {
        VariableUtils().SetFlag(rFlag, false, it.first->Elements());
    }

    // set entity flags for evaluated model parts
    for (auto& p_model_part : rEvaluatedModelParts) {
        VariableUtils().SetFlag(rFlag, true, p_model_part->Elements());
    }

    // get number of overlapping entities
    OptimizationUtils::SensitivityVariableModelPartsListMap reversed_map;
    OptimizationUtils::ReverseSensitivityModelPartVariablesListMap(reversed_map, rSensitivityModelPartVariableInfo);
    for (const auto& it : reversed_map) {

        IndexType number_of_common_entities = 0;
        for (const auto& p_model_part : it.second) {
            number_of_common_entities += OptimizationUtils::GetNumberOfContainerItemsWithFlag(p_model_part->Elements(), p_model_part->GetCommunicator().GetDataCommunicator(), rFlag);
        }

        std::visit([&](auto&& r_variable) {
                KRATOS_ERROR_IF(number_of_common_entities == 0)
                    << "No common entities between evaluated and sensitivity "
                       "model parts found for sensitivity variable "
                    << r_variable->Name() << ".\n";
        }, it.first);
    }

    // clear all the sensitivity variables for nodes. Here we assume there are
    // no overlapping regions in Elements and/or Conditions between provided rSensitivityModelParts hence, SetValue is
    // used in Elements and/or Condtions. Nodal sensitivities are added so that common nodes between two model parts
    // will have correct sensitivities.
    for (const auto& it_var_1 : rUsedNodalSensitivityVariables) {
        std::visit([&](auto&& p_var_1) {
            for (const auto& it : reversed_map) {
                std::visit([&](auto&& p_var_2) {
                    if (*p_var_1 == *p_var_2) {
                        for (const auto p_model_part : it.second) {
                            VariableUtils().SetNonHistoricalVariableToZero(*p_var_1, p_model_part->Nodes());
                        }
                    }
                }, it.first);
            }
        }, it_var_1);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos