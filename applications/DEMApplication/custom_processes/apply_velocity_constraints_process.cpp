
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//

// System includes

// External includes

// Project includes
#include "apply_velocity_constraints_process.hpp"
#include "utilities/parallel_utilities.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{
    /* Public functions *******************************************************/
    ApplyVelocityConstraintsProcess::ApplyVelocityConstraintsProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process(Flags()) , mrModelPart(rModelPart), mParameters(rParameters), mInterval(rParameters)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "help"                 : "This process applies constraints to the particles in a certain submodelpart, for a certain time interval",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "variable_name"        : "object",
                "interval"             : [0.0, 1e30],
                "constrained"          : [true,true,true],
                "value"                : [10.0, "3*t", "x+y"],
                "table"                : [0, 0, 0]
            } )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVelocityFunctions.clear();
        mpVelocityTable.clear();

        for(int i=0; i<3; i++) {
            mVelocityIsConstrained[i] = rParameters["constrained"][i].GetBool();
            if(rParameters["value"][i].IsNull()) {
                mVelocityValueIsNumeric[i] = true;
                mVelocityValues[i] = 0.0;
                mVelocityFunctions.push_back(GenericFunctionUtility("0.0"));
                // because I can't construct an array_1d of these
            } else {
                if(rParameters["value"][i].IsNumber()) {
                    mVelocityValueIsNumeric[i] = true;
                    mVelocityValues[i] = rParameters["value"][i].GetDouble();
                    mVelocityFunctions.push_back(GenericFunctionUtility("0.0"));
                    // because I can't construct an array_1d of these
                } else {
                    mVelocityValueIsNumeric[i] = false;
                    mVelocityFunctions.push_back(GenericFunctionUtility(rParameters["value"][i].GetString()));
                }
            }

            if(rParameters["table"][i].IsNull()) {
                mVelocityTableId[i] = 0;
            } else {
                mVelocityTableId[i] = rParameters["table"][i].GetInt();
            }
            mpVelocityTable.push_back(mrModelPart.pGetTable(mVelocityTableId[i]));
            // because I can't construct an array_1d of these
        }
        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ApplyVelocityConstraintsProcess::~ApplyVelocityConstraintsProcess() {}

    void ApplyVelocityConstraintsProcess::Execute() {}

    void ApplyVelocityConstraintsProcess::ExecuteInitialize()
    {
        KRATOS_TRY;

        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            array_1d<double, 3>& vel = rElement.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            for(int i=0; i<3; i++) {

                double velocity_value = 0.0;
                if(mVelocityValueIsNumeric[i]) {
                    velocity_value = mVelocityValues[i];
                }
                vel[i] = velocity_value;
            }
        });

        KRATOS_CATCH("");
    }

    void ApplyVelocityConstraintsProcess::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;

        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            array_1d<double, 3>& vel = rElement.GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

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

            for(int i=0; i<3; i++) {
                if (mVelocityTableId[i] != 0) {
                    vel[i] = mpVelocityTable[i]->GetValue(time);
                } else {
                    if(mVelocityIsConstrained[i]) {
                        double velocity_value = 0.0;
                        if(mVelocityValueIsNumeric[i]) {
                            velocity_value = mVelocityValues[i];
                        } else {
                            velocity_value = mVelocityFunctions[i].CallFunction(rElement.GetGeometry()[0].X(), rElement.GetGeometry()[0].Y(), rElement.GetGeometry()[0].Z(), time);
                        }
                        vel[i] = velocity_value;
                    }
                }
            }
        });

        KRATOS_CATCH("");
    }

    void ApplyVelocityConstraintsProcess::ExecuteFinalizeSolutionStep()
    {
        KRATOS_TRY;

        const double time = mrModelPart.GetProcessInfo()[TIME];
        if(!mInterval.IsInInterval(time)) return;

        block_for_each(mrModelPart.Elements(), [&](Element& rElement)
        {
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, false);
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, false);
            rElement.GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, false);
            rElement.GetGeometry()[0].pGetDof(VELOCITY_X)->FreeDof();
            rElement.GetGeometry()[0].pGetDof(VELOCITY_Y)->FreeDof();
            rElement.GetGeometry()[0].pGetDof(VELOCITY_Z)->FreeDof();
        });

        KRATOS_CATCH("");
    }

    std::string ApplyVelocityConstraintsProcess::Info() const
    {
        return "ApplyVelocityConstraintsProcess";
    }

    void ApplyVelocityConstraintsProcess::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyVelocityConstraintsProcess";
    }

    void ApplyVelocityConstraintsProcess::PrintData(std::ostream& rOStream) const {}

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ApplyVelocityConstraintsProcess& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }
}
