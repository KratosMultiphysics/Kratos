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
#include "containers/variable.h"
#include "containers/variable_data.h"
#include "includes/kratos_components.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(ModelPart& rModelPart,
                                                                       const Parameters& ThisParameters)
    : mrModelPart(rModelPart)
{
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("value"))
        << "Missing 'value' parameter in the parameters of '"
        << GeoApplyConstantScalarValueProcess::Info() << "'" << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("variable_name"))
        << "Missing 'variable_name' parameter in the parameters of '"
        << GeoApplyConstantScalarValueProcess::Info() << "'" << std::endl;

    mVariableName = ThisParameters["variable_name"].GetString();
    mIsFixed      = ThisParameters.Has("is_fixed") ? ThisParameters["is_fixed"].GetBool() : false;

    KRATOS_ERROR_IF(mIsFixed && !KratosComponents<Variable<double>>::Has(mVariableName))
        << "It is not possible to fix the variable '" << mVariableName
        << "' which is not of type Variable<double>.\n";

    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<VariableData>::Get(mVariableName)))
        << "Trying to fix a variable that is not in the rModelPart - variable name is "
        << mVariableName << std::endl;

    if (KratosComponents<Variable<double>>::Has(mVariableName)) {
        mDoubleValue = ThisParameters["value"].GetDouble();
    } else if (KratosComponents<Variable<int>>::Has(mVariableName)) {
        mIntValue = ThisParameters["value"].GetInt();
    } else if (KratosComponents<Variable<bool>>::Has(mVariableName)) {
        mBoolValue = ThisParameters["value"].GetBool();
    }
}

void GeoApplyConstantScalarValueProcess::ExecuteInitialize()
{
    // Since DoF are by definition double variables, fixing DoF is only relevant for variables of type double.
    if (mIsFixed && KratosComponents<Variable<double>>::Has(mVariableName)) {
        VariableUtils().ApplyFixity(KratosComponents<Variable<double>>::Get(mVariableName), true,
                                    mrModelPart.Nodes());
    }

    if (KratosComponents<Variable<double>>::Has(mVariableName)) {
        VariableUtils().SetVariable(KratosComponents<Variable<double>>::Get(mVariableName),
                                    mDoubleValue, mrModelPart.Nodes());
    } else if (KratosComponents<Variable<int>>::Has(mVariableName)) {
        VariableUtils().SetVariable(KratosComponents<Variable<int>>::Get(mVariableName), mIntValue,
                                    mrModelPart.Nodes());
    } else if (KratosComponents<Variable<bool>>::Has(mVariableName)) {
        VariableUtils().SetVariable(KratosComponents<Variable<bool>>::Get(mVariableName),
                                    mBoolValue, mrModelPart.Nodes());
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

} // namespace Kratos
