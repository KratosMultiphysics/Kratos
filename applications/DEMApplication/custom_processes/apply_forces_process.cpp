#include "apply_forces_process.hpp"
#include "utilities/parallel_utilities.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    ApplyForcesProcess::ApplyForcesProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : mrModelPart(rModelPart), mParameters(rParameters), mInterval(rParameters)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "help"                 : "This process applies loads over the particles and walls in a certain submodelpart, for a certain time interval",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "object",
                "constrained"          : [true,true,true],
                    "value"            : [10.0, "3*t", "x+y"],
                    "table"            : [0, 0, 0],
                "interval"             : [0.0, 1e30]
            } )" );


        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        mForceFunctions.clear();
        mpForceTable.clear();

        for(int i=0; i<3; i++) {
            if(rParameters["value"][i].IsNull()) {
                mForceValueIsNumeric[i] = true;
                mForceValues[i] = 0.0;
                mForceFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            } else {
                if(rParameters["value"][i].IsNumber()) {
                    mForceValueIsNumeric[i] = true;
                    mForceValues[i] = rParameters["value"][i].GetDouble();
                    mForceFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
                } else {
                    mForceValueIsNumeric[i] = false;
                    mForceFunctions.push_back(GenericFunctionUtility(rParameters["value"][i].GetString()));
                }
            }

            if(rParameters["table"][i].IsNull()) {
                mForceTableId[i] = 0;
            } else {
                mForceTableId[i] = rParameters["table"][i].GetInt();
            }
            mpForceTable.push_back(mrModelPart.pGetTable(mForceTableId[i])); // because I can't construct an array_1d of these
        }
        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ApplyForcesProcess::~ApplyForcesProcess() {}

    void ApplyForcesProcess::Execute() {}

    void ApplyForcesProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;

        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            array_1d<double, 3>& force = rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
            for(int i=0; i<3; i++) {
                if (mForceTableId[i] != 0) {
                    force[i] = mpForceTable[i]->GetValue(time);
                } else {
                    double force_value = 0.0;
                    if(mForceValueIsNumeric[i]) {
                        force_value = mForceValues[i];
                    } else {
                        force_value = mForceFunctions[i].CallFunction(rElement.GetGeometry()[0].X(), rElement.GetGeometry()[0].Y(), rElement.GetGeometry()[0].Z(), time);
                    }
                    force[i] = force_value;
                }
            }
        });

        KRATOS_CATCH("");
    }

    void ApplyForcesProcess::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;

        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = ZeroVector(3);
        });

        KRATOS_CATCH("");
    }

    std::string ApplyForcesProcess::Info() const
    {
        return "ApplyForcesProcess";
    }

    void ApplyForcesProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyForcesProcess";
    }

    void ApplyForcesProcess::PrintData(std::ostream& rOStream) const {}

    /* External functions *****************************************************/
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyForcesProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
