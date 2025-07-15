// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#include "geo_apply_constant_scalar_value_process.h"
#include "containers/model.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(Model& rModel, Parameters ThisParameters)
    : GeoApplyConstantScalarValueProcess(
          rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{
}

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(ModelPart& rModelPart, Parameters ThisParameters)
    : Process(Flags()), mrModelPart(rModelPart)
{
    KRATOS_TRY

    // Some values need to be mandatorily prescribed since no meaningful default value exist.
    // So that an error is thrown if they don't exist
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("value"))
        << "Missing 'value' parameter in ThisParameters" << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("variable_name"))
        << "Missing 'variable_name' parameter in ThisParameters" << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("model_part_name"))
        << "Missing 'model_part_name' parameter in ThisParameters" << std::endl;

    auto default_parameters = Parameters(R"(
    {
        "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
        "variable_name"   : "PLEASE_PRESCRIBE_VARIABLE_NAME",
        "is_fixed"        : false,
        "value"           : 1.0
    }  )");
    // Now validate against defaults -- this also ensures no type mismatch
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName = ThisParameters["variable_name"].GetString();
    mIsFixed      = ThisParameters["is_fixed"].GetBool();
    if (KratosComponents<Variable<double>>::Has(mVariableName)) { // case of double variable
        mDoubleValue = ThisParameters["value"].GetDouble();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<double>>::Get(mVariableName)))
            << "Trying to fix a variable that is not in the rModelPart - variable name is "
            << mVariableName << std::endl;
    } else if (KratosComponents<Variable<int>>::Has(mVariableName)) { // Case of int variable
        mIntValue = ThisParameters["value"].GetInt();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<int>>::Get(mVariableName)))
            << "Trying to fix a variable that is not in the rModelPart - variable name is "
            << mVariableName << std::endl;
        KRATOS_ERROR_IF(mIsFixed)
            << "Sorry it is not possible to fix variables of type Variable<int>. Only double "
               "variables or vector components can be fixed"
            << std::endl;
    } else if (KratosComponents<Variable<bool>>::Has(mVariableName)) { // Case of bool variable
        mBoolValue = ThisParameters["value"].GetBool();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<bool>>::Get(mVariableName)))
            << "Trying to fix a variable that is not in the rModelPart - variable name is "
            << mVariableName << std::endl;
        KRATOS_ERROR_IF(mIsFixed)
            << "Sorry it is not possible to fix variables of type Variable<bool>. Only double "
               "variables or vector components can be fixed"
            << std::endl;
    }

    KRATOS_CATCH("");
}

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(ModelPart& rModelPart,
                                                                       const Variable<double>& rVariable,
                                                                       const double DoubleValue,
                                                                       const bool   IsFixed)
    : mrModelPart(rModelPart), mDoubleValue(DoubleValue), mIsFixed(IsFixed)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(rVariable))
        << "Trying to fix a variable that is not in the rModelPart - variable name is " << rVariable
        << std::endl;

    mVariableName = rVariable.Name();

    KRATOS_CATCH("");
}

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(ModelPart& rModelPart,
                                                                       const Variable<int>& rVariable,
                                                                       const int IntValue,
                                                                       bool      IsFixed)
    : mrModelPart(rModelPart), mIntValue(IntValue), mIsFixed(IsFixed)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mIsFixed)
        << "Sorry it is not possible to fix variables of type Variable<int>. Only double variables "
           "or vector components can be fixed"
        << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(rVariable))
        << "Trying to fix a variable that is not in the rModelPart - variable name is " << rVariable
        << std::endl;

    mVariableName = rVariable.Name();

    KRATOS_CATCH("");
}

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(ModelPart& rModelPart,
                                                                       const Variable<bool>& rVariable,
                                                                       const bool BoolValue,
                                                                       bool       IsFixed)
    : mrModelPart(rModelPart), mBoolValue(BoolValue), mIsFixed(IsFixed)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mIsFixed)
        << "Sorry it is not possible to fix variables of type Variable<bool>. Only double "
           "variables or vector components can be fixed"
        << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(rVariable))
        << "Trying to fix a variable that is not in the rModelPart - variable name is " << rVariable
        << std::endl;

    mVariableName = rVariable.Name();

    KRATOS_CATCH("");
}

void GeoApplyConstantScalarValueProcess::ExecuteInitialize()
{
    KRATOS_TRY;

    if (KratosComponents<Variable<double>>::Has(mVariableName)) { // case of double variable
        InternalApplyValue<>(KratosComponents<Variable<double>>::Get(mVariableName), mIsFixed, mDoubleValue);
    } else if (KratosComponents<Variable<int>>::Has(mVariableName)) { // Case of int variable
        InternalApplyValueWithoutFixing<>(KratosComponents<Variable<int>>::Get(mVariableName), mIntValue);
    } else if (KratosComponents<Variable<bool>>::Has(mVariableName)) { // Case of bool variable
        InternalApplyValueWithoutFixing<>(KratosComponents<Variable<bool>>::Get(mVariableName), mBoolValue);
    } else {
        KRATOS_ERROR << "Not able to fix the variable. Attempting to fix variable: " << mVariableName
                     << std::endl;
    }

    KRATOS_CATCH("");
}

void GeoApplyConstantScalarValueProcess::ExecuteFinalize()
{
    // Since DoF are by definition double variables, freeing DoF is only relevant for variables of type double.
    if (mIsFixed && KratosComponents<Variable<double>>::Has(mVariableName)) {
        VariableUtils().ApplyFixity(KratosComponents<Variable<double>>::Get(mVariableName), false,
                                    mrModelPart.Nodes());
    }
}

template <class TVarType>
void GeoApplyConstantScalarValueProcess::InternalApplyValue(const TVarType&               rVariable,
                                                            const bool                    ToBeFixed,
                                                            const typename TVarType::Type Value)
{
    const std::size_t number_of_nodes = mrModelPart.Nodes().size();

    if (number_of_nodes != 0) {
        block_for_each(mrModelPart.Nodes(), [&](Node& rNode) {
            if constexpr (std::is_same<TVarType, Variable<double>>::value) { // For nodes
                if (ToBeFixed) {
                    rNode.Fix(rVariable);
                }
            }
            rNode.FastGetSolutionStepValue(rVariable) = Value;
        });
    }
}

template void GeoApplyConstantScalarValueProcess::InternalApplyValue<Variable<bool>>(
    const Variable<bool>& rVariable, const bool ToBeFixed, const bool Value);
template void GeoApplyConstantScalarValueProcess::InternalApplyValue<Variable<int>>(const Variable<int>& rVariable,
                                                                                    const bool ToBeFixed,
                                                                                    const int Value);
template void GeoApplyConstantScalarValueProcess::InternalApplyValue<Variable<double>>(
    const Variable<double>& rVariable, const bool ToBeFixed, const double Value);

template <class TVarType>
void GeoApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing(const TVarType& rVariable,
                                                                         const typename TVarType::Type Value)
{
    const std::size_t number_of_nodes = mrModelPart.Nodes().size();

    if (number_of_nodes != 0) {
        VariableUtils().SetVariable(rVariable, Value, mrModelPart.Nodes());
    }
}

template void GeoApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing<Variable<bool>>(
    const Variable<bool>& rVariable, const bool Value);
template void GeoApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing<Variable<int>>(
    const Variable<int>& rVariable, const int Value);
template void GeoApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing<Variable<double>>(
    const Variable<double>& rVariable, const double Value);

} // namespace Kratos
