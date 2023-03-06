//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <sstream>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes

// Include base h
#include "optimization_utils.h"

namespace Kratos
{

bool OptimizationUtils::IsVariableInList(
    const SensitivityFieldVariableTypes& rVariable,
    const std::vector<SensitivityFieldVariableTypes>& rVariablesList)
{
    bool is_found = false;
    std::visit([&](auto&& r_variable) {
        is_found = std::find_if(rVariablesList.begin(), rVariablesList.end(), [&](auto& rVariableListVariableVariant) {
            bool is_found = false;
            std::visit([&](auto&& r_variable_list_variable) {
                is_found = *(r_variable) == *(r_variable_list_variable);
            }, rVariableListVariableVariant);

            return is_found;
        }) != rVariablesList.end();
    }, rVariable);

    return is_found;
}

template<class TContainerType>
IndexType OptimizationUtils::GetNumberOfContainerItemsWithFlag(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator,
    const Flags& rFlag,
    const bool FlagValue)
{
    KRATOS_TRY

    return rDataCommunicator.SumAll(block_for_each<SumReduction<IndexType>>(rContainer, [&](const auto& rEntity) {
        return rEntity.Is(rFlag) == FlagValue;
    }));

    KRATOS_CATCH("");
}

template<class TContainerType>
GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(
    const TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    if (rContainer.size() > 0) {
        const auto first_geometry_type = rContainer.begin()->GetGeometry().GetGeometryType();
        const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity) -> IndexType {
            return first_geometry_type == rEntity.GetGeometry().GetGeometryType();
        });

        if (rDataCommunicator.AndReduceAll(local_value)) {
            return first_geometry_type;
        } else {
            return GeometryData::KratosGeometryType::Kratos_generic_type;
        }
    } else {
        return GeometryData::KratosGeometryType::Kratos_generic_type;
    }

    KRATOS_CATCH("")
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAllContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MinReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.AndReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(
    const TContainerType& rContainer,
    const Variable<TDataType>& rVariable,
    const DataCommunicator& rDataCommunicator)
{
    KRATOS_TRY

    const bool local_value = block_for_each<MaxReduction<bool>>(rContainer, [&](const auto& rEntity){
        return rEntity.GetProperties().Has(rVariable);
    });

    return rDataCommunicator.OrReduceAll(local_value);

    KRATOS_CATCH("");
}

template<class TContainerType>
void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(
    ModelPart& rModelPart,
    TContainerType& rContainer)
{
    KRATOS_TRY

    auto properties_id = block_for_each<MaxReduction<IndexType>>(rModelPart.PropertiesArray(), [](const auto pProperties) {
        return pProperties->Id();
    });

    // creation of properties is done in serial
    for (auto& r_entity : rContainer) {
        auto p_properties = rModelPart.CreateNewProperties(++properties_id);
        const auto& element_properties = r_entity.GetProperties();
        *p_properties = element_properties;
        p_properties->SetId(properties_id);
        r_entity.SetProperties(p_properties);
    }

    KRATOS_CATCH("");
}

template<>
IndexType OptimizationUtils::GetVariableDimension(
    const Variable<double>& rVariable,
    const IndexType DomainSize)
{
    return 1;
}

template<>
IndexType OptimizationUtils::GetVariableDimension(
    const Variable<array_1d<double, 3>>& rVariable,
    const IndexType DomainSize)
{
    return DomainSize;
}

void OptimizationUtils::CopySolutionStepVariablesList(
    ModelPart& rDestinationModelPart,
    const ModelPart& rOriginModelPart)
{
    rDestinationModelPart.GetNodalSolutionStepVariablesList() = rOriginModelPart.GetNodalSolutionStepVariablesList();
}

void OptimizationUtils::ActivateEntitiesAndCheckOverlappingRegions(
    const std::vector<ModelPart*>& rEvaluatedModelParts,
    const SensitivityModelPartVariablesListMap& rSensitivityModelPartVariableInfo,
    const Flags& rActivatedFlag,
    const std::vector<SensitivityFieldVariableTypes>& rAllowedNodalSensitivityVariables,
    const std::vector<SensitivityFieldVariableTypes>& rAllowedConditionSensitivityVariables,
    const std::vector<SensitivityFieldVariableTypes>& rAllowedElementSensitivityVariables)
{
    KRATOS_TRY

    const bool check_nodal_entities = rAllowedNodalSensitivityVariables.size() > 0;
    const bool check_condition_entities = rAllowedConditionSensitivityVariables.size() > 0;
    const bool check_element_entities = rAllowedElementSensitivityVariables.size() > 0;

    // set all the sensitivity model part rActivatedFlag flags to false to avoid counting previously
    // initialized rActivatedFlag flags to be counted as active.
    for (auto& it : rSensitivityModelPartVariableInfo) {
        auto& r_sensitivity_model_part = *(it.first);
        if (check_nodal_entities) {
            VariableUtils().SetFlag(rActivatedFlag, false, r_sensitivity_model_part.Nodes());
        }

        if (check_condition_entities) {
            VariableUtils().SetFlag(rActivatedFlag, false, r_sensitivity_model_part.Conditions());
        }

        if (check_element_entities) {
            VariableUtils().SetFlag(rActivatedFlag, false, r_sensitivity_model_part.Elements());
        }
    }

    // setting active flag for nodes and elements of the evaluated model parts
    for (auto p_model_part: rEvaluatedModelParts) {
        if (check_nodal_entities) {
            VariableUtils().SetFlag(rActivatedFlag, true, p_model_part->Nodes());
        }

        if (check_condition_entities) {
            VariableUtils().SetFlag(rActivatedFlag, true, p_model_part->Conditions());
        }

        if (check_element_entities) {
            VariableUtils().SetFlag(rActivatedFlag, true, p_model_part->Elements());
        }
    }

    // Calculating overlaping number of entities between evaluated model parts and sensitivity model parts.
    std::unordered_map<SensitivityFieldVariableTypes, IndexType> overlaping_number_of_entities_map;
    for (auto& it : rSensitivityModelPartVariableInfo) {
        const auto& r_sensitivity_model_part = *(it.first);
        const auto& r_data_communicator = r_sensitivity_model_part.GetCommunicator().GetDataCommunicator();
        for (auto& r_variable : it.second) {
            if (IsVariableInList(r_variable, rAllowedNodalSensitivityVariables)) {
                std::visit([&](auto&& r_variable) {
                    overlaping_number_of_entities_map[r_variable] += OptimizationUtils::GetNumberOfContainerItemsWithFlag(r_sensitivity_model_part.Nodes(), r_data_communicator, rActivatedFlag);
                }, r_variable);
            } else if (IsVariableInList(r_variable, rAllowedConditionSensitivityVariables)) {
                std::visit([&](auto&& r_variable) {
                    overlaping_number_of_entities_map[r_variable] += OptimizationUtils::GetNumberOfContainerItemsWithFlag(r_sensitivity_model_part.Conditions(), r_data_communicator, rActivatedFlag);
                }, r_variable);
            } else if (IsVariableInList(r_variable, rAllowedElementSensitivityVariables)) {
                std::visit([&](auto&& r_variable) {
                    overlaping_number_of_entities_map[r_variable] += OptimizationUtils::GetNumberOfContainerItemsWithFlag(r_sensitivity_model_part.Elements(), r_data_communicator, rActivatedFlag);
                }, r_variable);
            } else {
                    std::stringstream msg;
                    std::visit([&](auto&& r_variable) {
                        msg << "Unsupported sensitivity w.r.t. " << r_variable->Name()
                            << " requested for " << r_sensitivity_model_part.FullName() << "."
                            << "\nFollowings are supported sensitivity variables:";

                        if (check_nodal_entities) {
                            for (auto& r_supported_variable : rAllowedNodalSensitivityVariables) {
                                std::visit([&](auto&& r_variable) {
                                    msg << "\n\t" << r_variable->Name();
                                }, r_supported_variable);
                            }
                        }

                        if (check_condition_entities) {
                            for (auto& r_supported_variable : rAllowedConditionSensitivityVariables) {
                                std::visit([&](auto&& r_variable) {
                                    msg << "\n\t" << r_variable->Name();
                                }, r_supported_variable);
                            }
                        }

                        if (check_element_entities) {
                            for (auto& r_supported_variable : rAllowedElementSensitivityVariables) {
                                std::visit([&](auto&& r_variable) {
                                    msg << "\n\t" << r_variable->Name();
                                }, r_supported_variable);
                            }
                        }
                    }, r_variable);

                    KRATOS_ERROR << msg.str();
            }
        }
    }

    // Checking whether there is an overlap for every sensitivity variable
    for (const auto& it : overlaping_number_of_entities_map) {
        if (it.second == 0) {
            std::stringstream msg;

            std::visit([&](auto&& r_variable) {
                msg << "No overlapping entities found for " << r_variable->Name()
                << " between evaluated model parts and sensitivity model parts. Followings are the evaluated model part names:";
                for (auto p_model_part: rEvaluatedModelParts) {
                    msg << "\n\t" << p_model_part->FullName();
                }

                msg << "\nFollowings are the sensitivity model parts with requested sensitivity variables:";
                for (auto& it : rSensitivityModelPartVariableInfo) {
                    msg << "\n\t" << it.first->FullName() << ": ";
                    for (auto& r_variable : it.second) {
                        std::visit([&](auto&& r_variable) {
                            msg << "\n\t\t" << r_variable->Name();
                        }, r_variable);
                    }
                }
            }, it.first);

            KRATOS_ERROR << msg.str();
        }
    }

    KRATOS_CATCH("");
}

// template instantiations
template IndexType OptimizationUtils::GetNumberOfContainerItemsWithFlag(const ModelPart::NodesContainerType&, const DataCommunicator&, const Flags&, const bool);
template IndexType OptimizationUtils::GetNumberOfContainerItemsWithFlag(const ModelPart::ConditionsContainerType&, const DataCommunicator&, const Flags&, const bool);
template IndexType OptimizationUtils::GetNumberOfContainerItemsWithFlag(const ModelPart::ElementsContainerType&, const DataCommunicator&, const Flags&, const bool);

template GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ConditionsContainerType&, const DataCommunicator&);
template GeometryData::KratosGeometryType OptimizationUtils::GetContainerEntityGeometryType(const ModelPart::ElementsContainerType&, const DataCommunicator&);

template void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ConditionsContainerType&);
template void OptimizationUtils::CreateEntitySpecificPropertiesForContainer(ModelPart&, ModelPart::ElementsContainerType&);

template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAllContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<double>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ConditionsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);
template bool OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties(const ModelPart::ElementsContainerType&, const Variable<array_1d<double, 3>>&, const DataCommunicator& rDataCommunicator);

}