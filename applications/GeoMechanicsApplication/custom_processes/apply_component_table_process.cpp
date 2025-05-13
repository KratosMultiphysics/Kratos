// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "apply_component_table_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

#include <boost/algorithm/string.hpp>

namespace Kratos
{

ApplyComponentTableProcess::ApplyComponentTableProcess(ModelPart& rModelPart, Parameters ProcessSettings)
    : Process(Flags()), mrModelPart(rModelPart)
{
    KRATOS_TRY

    // only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "value" : 1.0,
                "table" : 1
            }  )");

    // Some values need to be mandatorily prescribed since no meaningful default value exist.
    // For this reason try accessing to them So that an error is thrown if they don't exist
    ProcessSettings["table"];
    ProcessSettings["variable_name"];
    ProcessSettings["model_part_name"];

    mIsFixedProvided = ProcessSettings.Has("is_fixed");
    // Now validate against defaults -- this also ensures no type mismatch
    ProcessSettings.ValidateAndAssignDefaults(default_parameters);

    mVariableName = ProcessSettings["variable_name"].GetString();
    mIsFixed      = ProcessSettings["is_fixed"].GetBool();
    mInitialValue = ProcessSettings["value"].GetDouble();

    unsigned int TableId = ProcessSettings["table"].GetInt();
    mpTable              = rModelPart.pGetTable(TableId);
    KRATOS_ERROR_IF_NOT(rModelPart.GetProcessInfo().Has(TIME_UNIT_CONVERTER))
        << "TIME_UNIT_CONVERTER not found in ProcessInfo for the ApplyComponentTableProcess" << std::endl;
    mTimeUnitConverter = rModelPart.GetProcessInfo()[TIME_UNIT_CONVERTER];

    KRATOS_CATCH("")
}

void ApplyComponentTableProcess::ExecuteInitialize()
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mVariableName))
        << "Variable '" << mVariableName << "' is not defined!" << std::endl;

    const auto& var = KratosComponents<Variable<double>>::Get(mVariableName);

    auto variable_name_1 = mpTable->NameOfX();
    boost::to_upper(variable_name_1);
    if (variable_name_1 == "TIME") {
        block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
            if (mIsFixed) rNode.Fix(var);
            else if (mIsFixedProvided) rNode.Free(var);
            rNode.FastGetSolutionStepValue(var) = mInitialValue;
        });
    } else if (variable_name_1 == "X") {
        block_for_each(mrModelPart.Nodes(), [&var, this](auto& rNode) {
            if (mIsFixed) rNode.Fix(var);
            else if (mIsFixedProvided) rNode.Free(var);
            rNode.FastGetSolutionStepValue(var) = mpTable->GetValue(rNode.X());
        });
    } else {
        KRATOS_ERROR
            << "Failed to initialize ApplyComponentTableProcess: got unknown table variable '"
            << variable_name_1 << "'";
    }

    KRATOS_CATCH("")
}

void ApplyComponentTableProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    const auto& var = KratosComponents<Variable<double>>::Get(mVariableName);

    auto variable_name_1 = mpTable->NameOfX();
    boost::to_upper(variable_name_1);
    if (variable_name_1 == "TIME") {
        const double Time  = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
        const double value = mpTable->GetValue(Time);
        block_for_each(mrModelPart.Nodes(),
                       [&var, &value](Node& rNode) { rNode.FastGetSolutionStepValue(var) = value; });
    }

    KRATOS_CATCH("")
}

void ApplyComponentTableProcess::ExecuteFinalize()
{
    if (mIsFixed) {
        VariableUtils().ApplyFixity(KratosComponents<Variable<double>>::Get(mVariableName), false,
                                    mrModelPart.Nodes());
    }
}

std::string ApplyComponentTableProcess::Info() const { return "ApplyComponentTableProcess"; }

} // namespace Kratos
