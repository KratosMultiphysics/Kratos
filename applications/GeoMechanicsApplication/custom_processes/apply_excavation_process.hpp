// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Lorenzo Gracia,
//                   Aron Noordam,
//                   Vahid Galavi
//
#pragma once

// Project includes
#include "includes/model_part.h"
#include "includes/element.h"
#include "utilities/variable_utils.h"

// Application includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"


namespace Kratos
{

class ApplyExcavationProcess : public Process
{
  public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyExcavationProcess);


    ApplyExcavationProcess(ModelPart& model_part,
                           Parameters Settings) : Process(Flags()), mrModelPart(model_part)
    {
        KRATOS_TRY
        mDeactivateSoilPart = Settings["deactivate_soil_part"].GetBool();
        KRATOS_CATCH("")
    }


    ~ApplyExcavationProcess() override = default;


    void ExecuteInitialize() override
    {
        KRATOS_TRY

        VariableUtils{}.SetFlag(ACTIVE, !mDeactivateSoilPart, mrModelPart.Elements());

        if (mDeactivateSoilPart) {
            block_for_each(mrModelPart.Elements(), [](Element& rElement) {
                rElement.ResetConstitutiveLaw();
            });
        } else {
            VariableUtils{}.SetFlag(ACTIVE, true, mrModelPart.Nodes());
        }

        VariableUtils{}.SetFlag(ACTIVE, !mDeactivateSoilPart, mrModelPart.Conditions());

        KRATOS_CATCH("")
    }


  private:
    ModelPart& mrModelPart;
    bool mDeactivateSoilPart;
};

}