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
//                   Richard Faasse
//

#include "geo_apply_constant_scalar_value_process.h"
#include "containers/variable.h"
#include "containers/variable_data.h"
#include "includes/kratos_components.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

#include <string>

namespace Kratos
{
using namespace std::string_literals;

GeoApplyConstantScalarValueProcess::GeoApplyConstantScalarValueProcess(ModelPart& rModelPart,
                                                                       const Parameters& rParameters)
    : mrModelPart(rModelPart)
{
    KRATOS_ERROR_IF_NOT(rParameters.Has("value"))
        << "Missing 'value' parameter in the parameters of '"
        << GeoApplyConstantScalarValueProcess::Info() << "'" << std::endl;
    KRATOS_ERROR_IF_NOT(rParameters.Has("variable_name"))
        << "Missing 'variable_name' parameter in the parameters of '"
        << GeoApplyConstantScalarValueProcess::Info() << "'" << std::endl;

    mVariableName = rParameters["variable_name"].GetString();
    if (rParameters.Has("is_fixed")) mIsFixed = rParameters["is_fixed"].GetBool();

    KRATOS_ERROR_IF(mIsFixed && !KratosComponents<Variable<double>>::Has(mVariableName))
        << "It is not possible to fix the variable '" << mVariableName
        << "' which is not of type Variable<double>.\n";

    KRATOS_ERROR_IF_NOT(rModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents<VariableData>::Get(mVariableName)))
        << "Trying to fix a variable that is not in ModelPart '" << rModelPart.Name()
        << "' - variable name is " << mVariableName << std::endl;

    if (KratosComponents<Variable<double>>::Has(mVariableName)) {
        mDoubleValue = rParameters["value"].GetDouble();
    } else if (KratosComponents<Variable<int>>::Has(mVariableName)) {
        mIntValue = rParameters["value"].GetInt();
    } else if (KratosComponents<Variable<bool>>::Has(mVariableName)) {
        mBoolValue = rParameters["value"].GetBool();
    }
}

void GeoApplyConstantScalarValueProcess::ExecuteInitialize()
{
    // Since DoF are by definition double variables, fixing DoF is only relevant for variables of type double.
    if (mIsFixed && KratosComponents<Variable<double>>::Has(mVariableName)) {
        VariableUtils().ApplyFixity(KratosComponents<Variable<double>>::Get(mVariableName), true,
                                    mrModelPart.Nodes());
    }
}

void GeoApplyConstantScalarValueProcess::ExecuteInitializeSolutionStep()
{
    if (mIsInitialized)
        return; // Constant value process, execute once here to ensure correct total and incremental D.O.F. values

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
    mIsInitialized = true;
}

void GeoApplyConstantScalarValueProcess::ExecuteFinalize()
{
    // Since DoF are by definition double variables, freeing DoF is only relevant for variables of type double.
    if (mIsFixed && KratosComponents<Variable<double>>::Has(mVariableName)) {
        VariableUtils().ApplyFixity(KratosComponents<Variable<double>>::Get(mVariableName), false,
                                    mrModelPart.Nodes());
    }
}

std::string GeoApplyConstantScalarValueProcess::Info() const
{
    return "GeoApplyConstantScalarValueProcess"s;
}

} // namespace Kratos
