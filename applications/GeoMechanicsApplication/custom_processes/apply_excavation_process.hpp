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
//


// System includes
#include <cmath>
#include <iostream>
#include<string>

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes

#if !defined(KRATOS_GEO_APPLY_EXCAVATION_PROCESS )
#define  KRATOS_GEO_APPLY_EXCAVATION_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyExcavationProcess : public Process
{
  public:

    typedef std::size_t IndexType;
    typedef Table<double, double> TableType;

    KRATOS_CLASS_POINTER_DEFINITION(ApplyExcavationProcess);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyExcavationProcess(ModelPart&  model_part,
                           Parameters& rParameters) : Process(Flags()), mrModelPart(model_part)
    {
        KRATOS_TRY

        mDeactivateSoilPart =  rParameters["deactivate_soil_part"].GetBool();
        mModelPartName      =  rParameters["model_part_name"].GetString();

        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyExcavationProcess() override{}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfElements() > 0) {
            if (mDeactivateSoilPart) {
                // Deactivation of the existing parts:
                block_for_each(mrModelPart.Elements(), [&](Element& rElement) {
                    rElement.Set(ACTIVE, false);
                    rElement.ResetConstitutiveLaw();
                });
            } else {
                // Activation of the existing parts:
                block_for_each(mrModelPart.Elements(), [&](Element& rElement) {
                    rElement.Set(ACTIVE, true);
                });

                // Same nodes for both computing model part
                if (mrModelPart.NumberOfNodes() > 0) {
                    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode) {
                        rNode.Set(ACTIVE, true);
                    });
                }
            }
        }

        // Conditions
        if (mrModelPart.NumberOfConditions() > 0) {
            if (mDeactivateSoilPart) {
                block_for_each(mrModelPart.Conditions(), [&](Condition& rCondition) {
                    rCondition.Set(ACTIVE, false);
                });
            } else {
                block_for_each(mrModelPart.Conditions(), [&](Condition& rCondition) {
                    rCondition.Set(ACTIVE, true);
                });
            }
        }

        KRATOS_CATCH("")
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables

    ModelPart& mrModelPart;
    bool mDeactivateSoilPart;
    std::string mModelPartName;

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_CONSTRUCTION_UTILITIES defined */