#include "apply_kinematic_constraints_to_walls_process.hpp"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    ApplyKinematicConstraintsToWallsProcess::ApplyKinematicConstraintsToWallsProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process(Flags()) , mrModelPart(rModelPart), mParameters(rParameters), mInterval(rParameters)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "help"                 : "This process applies constraints to the rigid walls in a certain submodelpart, for a certain time interval",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "velocity_constraints_settings" : {
                    "constrained"          : [true,true,true],
                    "value"                : [10.0, "3*t", "x+y"],
                    "table"                : [0, 0, 0]
                },
                "angular_velocity_constraints_settings" : {
                    "constrained"          : [true,true,true],
                    "value"                : [10.0, "3*t", "x+y"],
                    "table"                : [0, 0, 0]
                },
                "interval"             : [0.0, 1e30]
            } )" );



        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVelocityFunctions.clear();
        mAngularVelocityFunctions.clear();

        mpVelocityTable.clear();
        mpAngularVelocityTable.clear();

        for(int i=0; i<3; i++) {
            mVelocityIsConstrained[i] = rParameters["velocity_constraints_settings"]["constrained"][i].GetBool();
            mAngularVelocityIsConstrained[i] = rParameters["angular_velocity_constraints_settings"]["constrained"][i].GetBool();
            if(rParameters["velocity_constraints_settings"]["value"][i].IsNull()) {
                mVelocityValueIsNumeric[i] = true;
                mVelocityValues[i] = 0.0;
                mVelocityFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                if(rParameters["velocity_constraints_settings"]["value"][i].IsNumber()) {
                    mVelocityValueIsNumeric[i] = true;
                    mVelocityValues[i] = rParameters["velocity_constraints_settings"]["value"][i].GetDouble();
                    mVelocityFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
                }
                else {
                    mVelocityValueIsNumeric[i] = false;
                    mVelocityFunctions.push_back(GenericFunctionUtility(rParameters["velocity_constraints_settings"]["value"][i].GetString()));
                }
            }

            if(rParameters["velocity_constraints_settings"]["table"][i].IsNull()) {
                mVelocityTableId[i] = 0;
            }
            else {
                mVelocityTableId[i] = rParameters["velocity_constraints_settings"]["table"][i].GetInt();
            }
            mpVelocityTable.push_back(mrModelPart.pGetTable(mVelocityTableId[i])); // because I can't construct an array_1d of these

            if(rParameters["angular_velocity_constraints_settings"]["value"][i].IsNull()) {
                mAngularVelocityValueIsNumeric[i] = true;
                mAngularVelocityValues[i] = 0.0;
                mAngularVelocityFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                if(rParameters["angular_velocity_constraints_settings"]["value"][i].IsNumber()) {
                    mAngularVelocityValueIsNumeric[i] = true;
                    mAngularVelocityValues[i] = rParameters["angular_velocity_constraints_settings"]["value"][i].GetDouble();
                    mAngularVelocityFunctions.push_back(GenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
                }
                else {
                    mAngularVelocityValueIsNumeric[i] = false;
                    mAngularVelocityFunctions.push_back(GenericFunctionUtility(rParameters["angular_velocity_constraints_settings"]["value"][i].GetString()));
                }
            }
            if(rParameters["angular_velocity_constraints_settings"]["table"][i].IsNull()) {
                mAngularVelocityTableId[i] = 0;
            }
            else {
                mAngularVelocityTableId[i] = rParameters["angular_velocity_constraints_settings"]["table"][i].GetInt();
            }
            mpAngularVelocityTable.push_back(mrModelPart.pGetTable(mAngularVelocityTableId[i])); // because I can't construct an array_1d of these
        }

        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ApplyKinematicConstraintsToWallsProcess::~ApplyKinematicConstraintsToWallsProcess()
    {

    }

    void ApplyKinematicConstraintsToWallsProcess::Execute()
    {
    }

    void ApplyKinematicConstraintsToWallsProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {

            array_1d<double, 3>& vel = rElement.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& ang_vel = rElement.GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

            if(mVelocityIsConstrained[0]) {
                rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, true);
                rElement.GetGeometry()[0].pGetDof(VELOCITY_X)->FixDof();
            }
            if(mVelocityIsConstrained[1]) {
                rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, true);
                rElement.GetGeometry()[0].pGetDof(VELOCITY_Y)->FixDof();
            }
            if(mVelocityIsConstrained[2]) {
                rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, true);
                rElement.GetGeometry()[0].pGetDof(VELOCITY_Z)->FixDof();
            }
            if(mAngularVelocityIsConstrained[0]) {
                rElement.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X, true);
                rElement.GetGeometry()[0].pGetDof(ANGULAR_VELOCITY_X)->FixDof();
            }
            if(mAngularVelocityIsConstrained[1]) {
                rElement.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y, true);
                rElement.GetGeometry()[0].pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
            }
            if(mAngularVelocityIsConstrained[2]) {
                rElement.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z, true);
                rElement.GetGeometry()[0].pGetDof(ANGULAR_VELOCITY_Z)->FixDof();
            }

            for(int i=0; i<3; i++) {
                if (mVelocityTableId[i] != 0) {
                    vel[i] = mpVelocityTable[i]->GetValue(time);
                }
                else {
                    if(mVelocityIsConstrained[i]) {
                        double velocity_value = 0.0;
                        if(mVelocityValueIsNumeric[i]) {
                            velocity_value = mVelocityValues[i];
                        }
                        else {
                            velocity_value = mVelocityFunctions[i].CallFunction(rElement.GetGeometry()[0].X(), rElement.GetGeometry()[0].Y(), rElement.GetGeometry()[0].Z(), time);
                        }
                        vel[i] = velocity_value;
                    }
                }

                if (mAngularVelocityTableId[i] != 0) {
                    ang_vel[i] = mpAngularVelocityTable[i]->GetValue(time);
                }
                else {
                    if(mAngularVelocityIsConstrained[i]) {
                        double angular_velocity_value = 0.0;
                        if(mAngularVelocityValueIsNumeric[i]) {
                            angular_velocity_value = mAngularVelocityValues[i];
                        }
                        else {
                            angular_velocity_value = mAngularVelocityFunctions[i].CallFunction(rElement.GetGeometry()[0].X(), rElement.GetGeometry()[0].Y(), rElement.GetGeometry()[0].Z(), time);
                        }
                        ang_vel[i] = angular_velocity_value;
                    }
                }
            }
        });

        KRATOS_CATCH("");
    }

    void ApplyKinematicConstraintsToWallsProcess::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {

            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, false);
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, false);
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, false);
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X, false);
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y, false);
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z, false);
            rElement.GetGeometry()[0].pGetDof(VELOCITY_X)->FreeDof();
            rElement.GetGeometry()[0].pGetDof(VELOCITY_Y)->FreeDof();
            rElement.GetGeometry()[0].pGetDof(VELOCITY_Z)->FreeDof();
            rElement.GetGeometry()[0].pGetDof(ANGULAR_VELOCITY_X)->FreeDof();
            rElement.GetGeometry()[0].pGetDof(ANGULAR_VELOCITY_Y)->FreeDof();
            rElement.GetGeometry()[0].pGetDof(ANGULAR_VELOCITY_Z)->FreeDof();
        });

            KRATOS_CATCH("");
    }

    std::string ApplyKinematicConstraintsToWallsProcess::Info() const
    {
        return "ApplyKinematicConstraintsToWallsProcess";
    }

    void ApplyKinematicConstraintsToWallsProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyKinematicConstraintsToWallsProcess";
    }

    void ApplyKinematicConstraintsToWallsProcess::PrintData(std::ostream& rOStream) const
    {
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyKinematicConstraintsToWallsProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
