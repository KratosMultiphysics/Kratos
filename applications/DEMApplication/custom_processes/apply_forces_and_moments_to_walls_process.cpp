#include "apply_forces_and_moments_to_walls_process.hpp"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    ApplyForcesAndMomentsToWallsProcess::ApplyForcesAndMomentsToWallsProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : mrModelPart(rModelPart), mParameters(rParameters), mInterval(rParameters)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "help"                 : "This process applies loads over the rigid walls in a certain submodelpart, for a certain time interval",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "force_settings" : {
                    "value"            : [10.0, "3*t", "x+y"],
                    "table"            : [0, 0, 0]
                },
                "moment_settings" : {
                    "value"            : [10.0, "3*t", "x+y"],
                    "table"            : [0, 0, 0]
                },
                "interval"             : [0.0, 1e30]
            } )" );



        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mForceFunctions.clear();
        mMomentFunctions.clear();

        mpForceTable.clear();
        mpMomentTable.clear();

        for(int i=0; i<3; i++) {
            if(rParameters["force_settings"]["value"][i].IsNull()) {
                mForceValueIsNumeric[i] = true;
                mForceValues[i] = 0.0;
                mForceFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                if(rParameters["force_settings"]["value"][i].IsNumber()) {
                    mForceValueIsNumeric[i] = true;
                    mForceValues[i] = rParameters["force_settings"]["value"][i].GetDouble();
                    mForceFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
                }
                else {
                    mForceValueIsNumeric[i] = false;
                    mForceFunctions.push_back(GenericFunctionUtility(rParameters["force_settings"]["value"][i].GetString()));
                }
            }

            if(rParameters["force_settings"]["table"][i].IsNull()) {
                mForceTableId[i] = 0;
            }
            else {
                mForceTableId[i] = rParameters["force_settings"]["table"][i].GetInt();
            }
            mpForceTable.push_back(mrModelPart.pGetTable(mForceTableId[i])); // because I can't construct an array_1d of these

            if(rParameters["moment_settings"]["value"][i].IsNull()) {
                mMomentValueIsNumeric[i] = true;
                mMomentValues[i] = 0.0;
                mMomentFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                if(rParameters["moment_settings"]["value"][i].IsNumber()) {
                    mMomentValueIsNumeric[i] = true;
                    mMomentValues[i] = rParameters["moment_settings"]["value"][i].GetDouble();
                    mMomentFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
                }
                else {
                    mMomentValueIsNumeric[i] = false;
                    mMomentFunctions.push_back(GenericFunctionUtility(rParameters["moment_settings"]["value"][i].GetString()));
                }
            }

            if(rParameters["moment_settings"]["table"][i].IsNull()) {
                mMomentTableId[i] = 0;
            }
            else {
                mMomentTableId[i] = rParameters["moment_settings"]["table"][i].GetInt();
            }
            mpMomentTable.push_back(mrModelPart.pGetTable(mMomentTableId[i])); // because I can't construct an array_1d of these
        }

        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ApplyForcesAndMomentsToWallsProcess::~ApplyForcesAndMomentsToWallsProcess()
    {

    }

    void ApplyForcesAndMomentsToWallsProcess::Execute()
    {
    }

    void ApplyForcesAndMomentsToWallsProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {

            array_1d<double, 3>& force = rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
            array_1d<double, 3>& moment = rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);

            for(int i=0; i<3; i++) {
                if (mForceTableId[i] != 0) {
                    force[i] = mpForceTable[i]->GetValue(time);
                }
                else {
                    double force_value = 0.0;
                    if(mForceValueIsNumeric[i]) {
                        force_value = mForceValues[i];
                    }
                    else {
                        force_value = mForceFunctions[i].CallFunction(rElement.GetGeometry()[0].X(), rElement.GetGeometry()[0].Y(), rElement.GetGeometry()[0].Z(), time);
                    }
                    force[i] = force_value;
                }

                if (mMomentTableId[i] != 0) {
                    moment[i] = mpMomentTable[i]->GetValue(time);
                }
                else {
                    double moment_value = 0.0;
                    if(mMomentValueIsNumeric[i]) {
                        moment_value = mMomentValues[i];
                    }
                    else {
                        moment_value = mMomentFunctions[i].CallFunction(rElement.GetGeometry()[0].X(), rElement.GetGeometry()[0].Y(), rElement.GetGeometry()[0].Z(), time);
                    }
                    moment[i] = moment_value;
                }
            }
        });

        KRATOS_CATCH("");
    }

    void ApplyForcesAndMomentsToWallsProcess::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = ZeroVector(3);
            rElement.GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT) = ZeroVector(3);
        });

        KRATOS_CATCH("");
    }

    std::string ApplyForcesAndMomentsToWallsProcess::Info() const
    {
        return "ApplyForcesAndMomentsToWallsProcess";
    }

    void ApplyForcesAndMomentsToWallsProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyForcesAndMomentsToWallsProcess";
    }

    void ApplyForcesAndMomentsToWallsProcess::PrintData(std::ostream& rOStream) const
    {
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyForcesAndMomentsToWallsProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
