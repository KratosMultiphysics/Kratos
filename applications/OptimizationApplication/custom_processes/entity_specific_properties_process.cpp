//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "entity_specific_properties_process.h"

namespace Kratos {

const Parameters EntitySpecificPropertiesProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
            "variables_list" : ["PLEASE_SPECIFY_VARIABLE_NAMES_LIST"],
            "echo_level"     : 0
        })");

    return default_parameters;
}

EntitySpecificPropertiesProcess::EntitySpecificPropertiesProcess(
    Model& rModel,
    Parameters rParameters)
    : mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mModelPartName = rParameters["model_part_name"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();

    for (const auto& r_variable_name : rParameters["variables_list"].GetStringArray()) {
        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            mPropertiesSpecificVariablePointersList.push_back(&(KratosComponents<Variable<double>>::Get(r_variable_name)));
        } else {
            KRATOS_ERROR << r_variable_name << " is not of type double variable.";
        }
    }

    KRATOS_CATCH("");
}

int EntitySpecificPropertiesProcess::Check()
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("");
}

void EntitySpecificPropertiesProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    if (!mIsPropertiesReplaced) {
        CreateEntitySpecificProperties(r_model_part.Conditions());
        CreateEntitySpecificProperties(r_model_part.Elements());
        mIsPropertiesReplaced = true;
    }

    UpdateEntitySpecificProperties(r_model_part.Conditions());
    UpdateEntitySpecificProperties(r_model_part.Elements());

    KRATOS_CATCH("");
}

std::string EntitySpecificPropertiesProcess::Info() const
{
    return std::string("EntitySpecificPropertiesProcess");
}

void EntitySpecificPropertiesProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void EntitySpecificPropertiesProcess::PrintData(std::ostream& rOStream) const
{
}

template<class ContainerType>
void EntitySpecificPropertiesProcess::CreateEntitySpecificProperties(ContainerType& rContainer)
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);

    // creation of properties is done in serial
    int properties_id = r_model_part.NumberOfProperties() + 1;
    for (auto& r_entity : rContainer) {
        auto p_properties = r_model_part.CreateNewProperties(properties_id++);
        const auto& element_properties = r_entity.GetProperties();
        *p_properties = element_properties;
        r_entity.SetProperties(p_properties);
    }

    // now make the variable space in entities
    block_for_each(rContainer, [&](auto& rEntity) {
        for (const auto p_var : mPropertiesSpecificVariablePointersList) {
            const auto& r_properties = rEntity.GetProperties();

            KRATOS_ERROR_IF_NOT(r_properties.Has(*p_var))
                << p_var->Name() << " is not found in the properties with id "
                << r_properties.Id() << " in model part " << r_model_part.FullName() << ".\n";

            rEntity.SetValue(*p_var, r_properties[*p_var]);
        }
    });

    KRATOS_CATCH("");
}

template<class ContainerType>
void EntitySpecificPropertiesProcess::UpdateEntitySpecificProperties(ContainerType& rContainer)
{
    KRATOS_TRY

    block_for_each(rContainer, [&](auto& rEntity) {
        for (const auto p_var : mPropertiesSpecificVariablePointersList) {
            rEntity.GetProperties()[*p_var] = rEntity.GetValue(*p_var);
        }
    });

    KRATOS_CATCH("");
}

} // namespace Kratos.
