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

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyConstantHydrostaticPressureProcess : public Process
{

public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantHydrostaticPressureProcess);

    ApplyConstantHydrostaticPressureProcess(ModelPart &model_part,
                                            Parameters rParameters) : Process(Flags()), mrModelPart(model_part)
    {
        KRATOS_TRY

        // only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
        {
            "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
            "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
            "is_fixed": false,
            "is_seepage" : false,
            "gravity_direction" : 2,
            "reference_coordinate" : 0.0,
            "specific_weight" : 10000.0,
            "pressure_tension_cut_off" : 0.0,
            "table" : 1
        }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["reference_coordinate"];
        rParameters["variable_name"];
        rParameters["model_part_name"];
        mIsFixedProvided = rParameters.Has("is_fixed");

        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mModelPartName       = rParameters["model_part_name"].GetString();
        mVariableName        = rParameters["variable_name"].GetString();
        mIsFixed             = rParameters["is_fixed"].GetBool();
        mIsSeepage           = rParameters["is_seepage"].GetBool();
        mGravityDirection    = rParameters["gravity_direction"].GetInt();
        mReferenceCoordinate = rParameters["reference_coordinate"].GetDouble();
        mSpecificWeight      = rParameters["specific_weight"].GetDouble();
        mPressureTensionCutOff = rParameters["pressure_tension_cut_off"].GetDouble();

        KRATOS_CATCH("")
    }

    ApplyConstantHydrostaticPressureProcess &operator=(const ApplyConstantHydrostaticPressureProcess&) = delete;
    ~ApplyConstantHydrostaticPressureProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        const Variable<double> &var = KratosComponents<Variable<double>>::Get(mVariableName);
        block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
            const double pressure = - PORE_PRESSURE_SIGN_FACTOR * mSpecificWeight * (mReferenceCoordinate - rNode.Coordinates()[mGravityDirection]);

            if (mIsSeepage) {
                if (pressure < PORE_PRESSURE_SIGN_FACTOR * mPressureTensionCutOff) {  // Before 0. was used i.s.o. the tension cut off value -> no effect in any test.
                    rNode.FastGetSolutionStepValue(var) = pressure;
                    if (mIsFixed) rNode.Fix(var);
                } else {
                    if (mIsFixedProvided) rNode.Free(var);
                }
            } else {
                rNode.FastGetSolutionStepValue(var) = std::min(pressure, PORE_PRESSURE_SIGN_FACTOR * mPressureTensionCutOff);
                if (mIsFixed) rNode.Fix(var);
                else if (mIsFixedProvided) rNode.Free(var);
            }
        });

        KRATOS_CATCH("")
    }

    const std::string& GetName() const
    {
        return mModelPartName;
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyConstantHydrostaticPressureProcess";
    }

protected:
    /// Member Variables
    ModelPart &mrModelPart;
    std::string mModelPartName;
    std::string mVariableName;
    bool mIsFixed;
    bool mIsFixedProvided;
    bool mIsSeepage;
    unsigned int mGravityDirection;
    double mReferenceCoordinate;
    double mSpecificWeight;
    double mPressureTensionCutOff;

};

}