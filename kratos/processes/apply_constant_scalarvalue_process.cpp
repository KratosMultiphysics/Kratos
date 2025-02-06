//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/variable_utils.h"
#include "processes/apply_constant_scalarvalue_process.h"

namespace Kratos
{
KRATOS_CREATE_LOCAL_FLAG(ApplyConstantScalarValueProcess, VARIABLE_IS_FIXED, 0)

/***********************************************************************************/
/***********************************************************************************/

ApplyConstantScalarValueProcess::ApplyConstantScalarValueProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : ApplyConstantScalarValueProcess(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()), ThisParameters)
{

}

/***********************************************************************************/
/***********************************************************************************/

ApplyConstantScalarValueProcess::ApplyConstantScalarValueProcess(
    ModelPart& rModelPart,
    Parameters ThisParameters
    ) : Process(Flags()) ,
        mrModelPart(rModelPart)
{
    KRATOS_TRY

    // Some values need to be mandatorily prescribed since no meaningful default value exist.
    // So that an error is thrown if they don't exist
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("value")) << "Missing 'value' parameter in ThisParameters" << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("variable_name")) << "Missing 'variable_name' parameter in ThisParameters" << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("model_part_name")) << "Missing 'model_part_name' parameter in ThisParameters" << std::endl;

    // Now validate against defaults -- this also ensures no type mismatch
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mMeshId = ThisParameters["mesh_id"].GetInt();
    mVariableName = ThisParameters["variable_name"].GetString();
    this->Set( VARIABLE_IS_FIXED, ThisParameters["is_fixed"].GetBool());

    if( KratosComponents<Variable<double>>::Has( mVariableName ) ) { //case of double variable
        mDoubleValue = ThisParameters["value"].GetDouble();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<double>>::Get(mVariableName))) << "Trying to fix a variable that is not in the rModelPart - variable name is " << mVariableName << std::endl;
    } else if( KratosComponents<Variable<int>>::Has( mVariableName ) ) { // Case of int variable
        mIntValue = ThisParameters["value"].GetInt();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<int>>::Get(mVariableName))) << "Trying to fix a variable that is not in the rModelPart - variable name is " << mVariableName << std::endl;
        KRATOS_ERROR_IF(this->Is(VARIABLE_IS_FIXED)) << "Sorry it is not possible to fix variables of type Variable<int>. Only double variables or vector components can be fixed" << std::endl;
    } else if( KratosComponents< Variable<bool> >::Has( mVariableName ) ) { // Case of bool variable
        mBoolValue = ThisParameters["value"].GetBool();
        KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<bool>>::Get(mVariableName))) << "Trying to fix a variable that is not in the rModelPart - variable name is " << mVariableName << std::endl;
        KRATOS_ERROR_IF(this->Is(VARIABLE_IS_FIXED)) << "Sorry it is not possible to fix variables of type Variable<bool>. Only double variables or vector components can be fixed" << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

ApplyConstantScalarValueProcess::ApplyConstantScalarValueProcess(
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const double DoubleValue,
    std::size_t MeshId,
    const Flags Options
    ) : Process(Options) ,
        mrModelPart(rModelPart),
        mDoubleValue(DoubleValue),
        mMeshId(MeshId)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(this->IsDefined(VARIABLE_IS_FIXED) == false) << "Please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)" << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(rVariable)) << "Trying to fix a variable that is not in the rModelPart - variable name is " << rVariable << std::endl;

    mVariableName = rVariable.Name();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

ApplyConstantScalarValueProcess::ApplyConstantScalarValueProcess(
    ModelPart& rModelPart,
    const Variable<int>& rVariable,
    const int IntValue,
    std::size_t MeshId,
    const Flags Options
    ) : Process(Options),
        mrModelPart(rModelPart),
        mIntValue(IntValue),
        mMeshId(MeshId)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(this->IsDefined(VARIABLE_IS_FIXED)) << "Please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)" << std::endl;
    KRATOS_ERROR_IF(this->Is(VARIABLE_IS_FIXED)) << "Sorry it is not possible to fix variables of type Variable<int>. Only double variables or vector components can be fixed" << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(rVariable)) << "Trying to fix a variable that is not in the rModelPart - variable name is " << rVariable << std::endl;

    mVariableName = rVariable.Name();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

ApplyConstantScalarValueProcess::ApplyConstantScalarValueProcess(
    ModelPart& rModelPart,
    const Variable< bool >& rVariable,
    const bool BoolValue,
    std::size_t MeshId,
    const Flags options
        ) : Process(options) ,
        mrModelPart(rModelPart),
        mBoolValue(BoolValue),
        mMeshId(MeshId)
{
    KRATOS_TRY;

    KRATOS_ERROR_IF_NOT(this->IsDefined(VARIABLE_IS_FIXED)) << "Please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)" << std::endl;
    KRATOS_ERROR_IF(this->Is(VARIABLE_IS_FIXED)) << "Sorry it is not possible to fix variables of type Variable<bool>. Only double variables or vector components can be fixed" << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(rVariable)) << "Trying to fix a variable that is not in the rModelPart - variable name is " << rVariable << std::endl;

    mVariableName = rVariable.Name();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ApplyConstantScalarValueProcess::ExecuteInitialize()
{
    KRATOS_TRY;
    const bool is_fixed = this->Is(VARIABLE_IS_FIXED);

    if( KratosComponents<Variable<double>>::Has( mVariableName ) ) { //case of double variable
        InternalApplyValue<>(KratosComponents< Variable<double> >::Get(mVariableName) , is_fixed, mDoubleValue);
    } else if( KratosComponents<Variable<int>>::Has( mVariableName ) ) { // Case of int variable
        InternalApplyValueWithoutFixing<>(KratosComponents<Variable<int>>::Get(mVariableName) , mIntValue);
    } else if( KratosComponents< Variable<bool> >::Has( mVariableName ) )  { // Case of bool variable
        InternalApplyValueWithoutFixing<>(KratosComponents<Variable<bool>>::Get(mVariableName), mBoolValue);
    } else {
        KRATOS_ERROR << "Not able to fix the variable. Attempting to fix variable: " << mVariableName << std::endl;
    }

    KRATOS_CATCH("");
}

void ApplyConstantScalarValueProcess::ExecuteFinalize() {
    if (this->Is(VARIABLE_IS_FIXED) && KratosComponents<Variable<double>>::Has(mVariableName)) {
        VariableUtils().ApplyFixity(
            KratosComponents<Variable<double>>::Get(mVariableName), false, mrModelPart.Nodes());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TVarType>
void ApplyConstantScalarValueProcess::InternalApplyValue(
    const TVarType& rVariable,
    const bool ToBeFixed,
    const typename TVarType::Type Value
    )
{
    const std::size_t number_of_nodes = mrModelPart.GetMesh(mMeshId).Nodes().size();

    if(number_of_nodes != 0) {
        block_for_each(mrModelPart.GetMesh(mMeshId).Nodes(), [&](Node& rNode){
            if constexpr (std::is_same<TVarType, Variable<double>>::value) { // For nodes
                if(ToBeFixed) {
                    rNode.Fix(rVariable);
                }
            }
            rNode.FastGetSolutionStepValue(rVariable) = Value;
        });
    }
}

template void ApplyConstantScalarValueProcess::InternalApplyValue<Variable<bool>>(
    const Variable<bool>& rVariable,
    const bool ToBeFixed,
    const bool Value
    );
template void ApplyConstantScalarValueProcess::InternalApplyValue<Variable<int>>(
    const Variable<int>& rVariable,
    const bool ToBeFixed,
    const int Value
    );
template void ApplyConstantScalarValueProcess::InternalApplyValue<Variable<double>>(
    const Variable<double>& rVariable,
    const bool ToBeFixed,
    const double Value
    );


/***********************************************************************************/
/***********************************************************************************/

template<class TVarType>
void ApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing(
    const TVarType& rVariable,
    const typename TVarType::Type Value
    )
{
    const std::size_t number_of_nodes = mrModelPart.GetMesh(mMeshId).Nodes().size();

    if(number_of_nodes != 0) {
        VariableUtils().SetVariable(rVariable, Value, mrModelPart.GetMesh(mMeshId).Nodes());
    }
}

template void ApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing<Variable<bool>>(
    const Variable<bool>& rVariable,
    const bool Value
    );
template void ApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing<Variable<int>>(
    const Variable<int>& rVariable,
    const int Value
    );
template void ApplyConstantScalarValueProcess::InternalApplyValueWithoutFixing<Variable<double>>(
    const Variable<double>& rVariable,
    const double Value
    );

/***********************************************************************************/
/***********************************************************************************/

const Parameters ApplyConstantScalarValueProcess::GetDefaultParameters() const
{
    return Parameters( R"(
    {
        "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
        "mesh_id"         : 0,
        "variable_name"   : "PLEASE_PRESCRIBE_VARIABLE_NAME",
        "is_fixed"        : false,
        "value"           : 1.0
    }  )" );
}

}
