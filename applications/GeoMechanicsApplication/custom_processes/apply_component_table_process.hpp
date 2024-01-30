// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

#include <boost/algorithm/string.hpp>

namespace Kratos
{

class ApplyComponentTableProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyComponentTableProcess);
    
    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double,double>;
    
    ApplyComponentTableProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "value" : 1.0,
                "table" : 1
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["table"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        mIsFixedProvided = rParameters.Has("is_fixed");
        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mIsFixed      = rParameters["is_fixed"].GetBool();
        mInitialValue = rParameters["value"].GetDouble();

        unsigned int TableId = rParameters["table"].GetInt();
        mpTable              = model_part.pGetTable(TableId);
        mTimeUnitConverter   = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
        
        KRATOS_CATCH("")
    }
    ApplyComponentTableProcess(const ApplyComponentTableProcess&) = delete;
    ApplyComponentTableProcess& operator=(const ApplyComponentTableProcess&) = delete;
    ~ApplyComponentTableProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
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
        }
        else if (variable_name_1 == "X") {
            block_for_each(mrModelPart.Nodes(), [&var, this](auto& node){
                node.FastGetSolutionStepValue(var) = mpTable->GetValue(node.X());
            });
        }
        else {
            KRATOS_ERROR << "Failed to initialize ApplyComponentTableProcess: got unknown table variable '"
                         << variable_name_1 << "'";
        }

        KRATOS_CATCH("")
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        const auto& var = KratosComponents<Variable<double>>::Get(mVariableName);

        auto variable_name_1 = mpTable->NameOfX();
        boost::to_upper(variable_name_1);
        if (variable_name_1 == "TIME") {
            const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
            const double value = mpTable->GetValue(Time);
            block_for_each(mrModelPart.Nodes(), [&var, &value](Node& rNode) {
                rNode.FastGetSolutionStepValue(var) = value;
            });
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyComponentTableProcess";
    }

private:
    /// Member Variables
    ModelPart& mrModelPart;
    std::string mVariableName;
    bool mIsFixed;
    bool mIsFixedProvided;
    double mInitialValue;
    TableType::Pointer mpTable;
    double mTimeUnitConverter;
};

}