#include "apply_moments_process.hpp"
#include "utilities/parallel_utilities.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    ApplyMomentsProcess::ApplyMomentsProcess(
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
        mMomentFunctions.clear();
        mpMomentTable.clear();

        for(int i=0; i<3; i++) {

            if(rParameters["value"][i].IsNull()) {
                mMomentValueIsNumeric[i] = true;
                mMomentValues[i] = 0.0;
                mMomentFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            } else {
                if(rParameters["value"][i].IsNumber()) {
                    mMomentValueIsNumeric[i] = true;
                    mMomentValues[i] = rParameters["value"][i].GetDouble();
                    mMomentFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
                } else {
                    mMomentValueIsNumeric[i] = false;
                    mMomentFunctions.push_back(GenericFunctionUtility(rParameters["value"][i].GetString()));
                }
            }

            if(rParameters["table"][i].IsNull()) {
                mMomentTableId[i] = 0;
            } else {
                mMomentTableId[i] = rParameters["table"][i].GetInt();
            }
            mpMomentTable.push_back(mrModelPart.pGetTable(mMomentTableId[i])); // because I can't construct an array_1d of these
        }
        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ApplyMomentsProcess::~ApplyMomentsProcess() {}

    void ApplyMomentsProcess::Execute() {}

    void ApplyMomentsProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;

        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            array_1d<double, 3>& moment = rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);

            for(int i=0; i<3; i++) {
                if (mMomentTableId[i] != 0) {
                    moment[i] = mpMomentTable[i]->GetValue(time);
                } else {
                    double moment_value = 0.0;
                    if(mMomentValueIsNumeric[i]) {
                        moment_value = mMomentValues[i];
                    } else {
                        moment_value = mMomentFunctions[i].CallFunction(rElement.GetGeometry()[0].X(), rElement.GetGeometry()[0].Y(), rElement.GetGeometry()[0].Z(), time);
                    }
                    moment[i] = moment_value;
                }
            }
        });

        KRATOS_CATCH("");
    }

    void ApplyMomentsProcess::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;

        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT) = ZeroVector(3);
        });

        KRATOS_CATCH("");
    }

    std::string ApplyMomentsProcess::Info() const
    {
        return "ApplyMomentsProcess";
    }

    void ApplyMomentsProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyMomentsProcess";
    }

    void ApplyMomentsProcess::PrintData(std::ostream& rOStream) const {}

    /* External functions *****************************************************/
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyMomentsProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
