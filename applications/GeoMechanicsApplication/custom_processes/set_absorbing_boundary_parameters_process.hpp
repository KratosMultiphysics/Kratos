// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class SetAbsorbingBoundaryParametersProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SetAbsorbingBoundaryParametersProcess);

    SetAbsorbingBoundaryParametersProcess(ModelPart& model_part, Parameters rParameters)
        : Process(Flags()), mrModelPart(model_part)
    {
        KRATOS_TRY

        // only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "absorbing_factors": [1.0,1.0],
                "virtual_thickness": 1e10,
                "skip_internal_forces": false
            }  )");

        // Some values need to be mandatory prescribed since no meaningful default value exist. For
        // this reason try accessing to them So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        // get absorbing factors
        mAbsorbingFactors.resize(2, false);
        mAbsorbingFactors(0) = rParameters["absorbing_factors"][0].GetDouble();
        mAbsorbingFactors(1) = rParameters["absorbing_factors"][1].GetDouble();

        // get virtual thickness
        mVirtualThickness = rParameters["virtual_thickness"].GetDouble();

        mSkipInternalForces = rParameters["skip_internal_forces"].GetBool();
        std::cout << "mSkipInternalForces: " << mSkipInternalForces << std::endl;

        KRATOS_CATCH("")
    }

    SetAbsorbingBoundaryParametersProcess(const SetAbsorbingBoundaryParametersProcess&) = delete;
    SetAbsorbingBoundaryParametersProcess& operator=(const SetAbsorbingBoundaryParametersProcess&) = delete;
    ~SetAbsorbingBoundaryParametersProcess() override = default;

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        block_for_each(mrModelPart.Conditions(), [&](Condition& rCondition) {
            rCondition.SetValue(ABSORBING_FACTORS, mAbsorbingFactors);
            rCondition.SetValue(VIRTUAL_THICKNESS, mVirtualThickness);
            rCondition.SetValue(SKIP_INTERNAL_FORCES, mSkipInternalForces);
        });

        KRATOS_CATCH("")
    }

    std::string Info() const override { return "SetAbsorbingBoundaryParametersProcess"; }

private:
    /// Member Variables
    ModelPart& mrModelPart;
    Vector     mAbsorbingFactors;
    double     mVirtualThickness;
	bool       mSkipInternalForces;
};

} // namespace Kratos