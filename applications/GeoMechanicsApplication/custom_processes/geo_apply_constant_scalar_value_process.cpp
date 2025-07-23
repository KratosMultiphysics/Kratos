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
#include "utilities/variable_utils.h"

namespace Kratos
{

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(ModelPart& rModelPart, Parameters ThisParameters)
    : mrModelPart(rModelPart)
{
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("value"))
        << "Missing 'value' parameter in ThisParameters" << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("variable_name"))
        << "Missing 'variable_name' parameter in ThisParameters" << std::endl;

    mVariableName = ThisParameters["variable_name"].GetString();
    mIsFixed      = ThisParameters.Has("is_fixed") ? ThisParameters["is_fixed"].GetBool() : false;
    if (KratosComponents<Variable<double>>::Has(mVariableName)) {
        mDoubleValue = ThisParameters["value"].GetDouble();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<double>>::Get(mVariableName)))
            << "Trying to fix a variable that is not in the rModelPart - variable name is "
            << mVariableName << std::endl;
    } else if (KratosComponents<Variable<int>>::Has(mVariableName)) {
        mIntValue = ThisParameters["value"].GetInt();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<int>>::Get(mVariableName)))
            << "Trying to fix a variable that is not in the rModelPart - variable name is "
            << mVariableName << std::endl;
        KRATOS_ERROR_IF(mIsFixed)
            << "Sorry it is not possible to fix variables of type Variable<int>. Only double "
               "variables or vector components can be fixed"
            << std::endl;
    } else if (KratosComponents<Variable<bool>>::Has(mVariableName)) {
        mBoolValue = ThisParameters["value"].GetBool();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<bool>>::Get(mVariableName)))
            << "Trying to fix a variable that is not in the rModelPart - variable name is "
            << mVariableName << std::endl;
        KRATOS_ERROR_IF(mIsFixed)
            << "Sorry it is not possible to fix variables of type Variable<bool>. Only double "
               "variables or vector components can be fixed"
            << std::endl;
    }
}

void GeoApplyConstantScalarValueProcess::ExecuteInitialize()
{
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
    if (!mrModelPart.Nodes().empty()) {
        if constexpr (std::is_same_v<TVarType, Variable<double>>) {
            if (ToBeFixed) {
                block_for_each(mrModelPart.Nodes(),
                               [&rVariable](Node& rNode) { rNode.Fix(rVariable); });
            }
        }

        block_for_each(mrModelPart.Nodes(), [&rVariable, &Value](Node& rNode) {
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
