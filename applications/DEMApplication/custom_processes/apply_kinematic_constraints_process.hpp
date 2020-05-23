//
//   Author:   Miguel Angel Celigueta
//

#if !defined(KRATOS_APPLY_KINEMATIC_CONSTRAINTS_PROCESS )
#define  KRATOS_APPLY_KINEMATIC_CONSTRAINTS_PROCESS

#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/interval_utility.h"
#include "utilities/python_function_callback_utility.h"
#include "DEM_application_variables.h"

namespace Kratos
{

class ApplyKinematicConstraintsProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyKinematicConstraintsProcess);


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyKinematicConstraintsProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mrModelPart(model_part), mParameters(rParameters), mInterval(rParameters)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "help"                 : "This process applies constraints to the particles in a certain submodelpart, for a certain time interval",
                "mesh_id"              : 0,
                "model_part_name"      : "please_specify_model_part_name",
                "velocity_constraints_settings" : {
                    "constrained"          : [true,true,true],
                    "value"                : [10.0, "3*t", "x+y"]
                },
                "angular_velocity_constraints_settings" : {
                    "constrained"          : [true,true,true],
                    "value"                : [10.0, "3*t", "x+y"]
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

        for(int i=0; i<3; i++) {
            mVelocityIsConstrained[i] = rParameters["velocity_constraints_settings"]["constrained"][i].GetBool();
            mAngularVelocityIsConstrained[i] = rParameters["angular_velocity_constraints_settings"]["constrained"][i].GetBool();
            if(rParameters["velocity_constraints_settings"]["value"][i].IsNumber()) {
                mVelocityValueIsNumeric[i] = true;
                mVelocityValues[i] = rParameters["velocity_constraints_settings"]["value"][i].GetDouble();
                mVelocityFunctions.push_back(PythonGenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                mVelocityValueIsNumeric[i] = false;
                mVelocityFunctions.push_back(PythonGenericFunctionUtility(rParameters["velocity_constraints_settings"]["value"][i].GetString()));
            }
            if(rParameters["angular_velocity_constraints_settings"]["value"][i].IsNumber()) {
                mAngularVelocityValueIsNumeric[i] = true;
                mAngularVelocityValues[i] = rParameters["angular_velocity_constraints_settings"]["value"][i].GetDouble();
                mAngularVelocityFunctions.push_back(PythonGenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                mAngularVelocityValueIsNumeric[i] = false;
                mAngularVelocityFunctions.push_back(PythonGenericFunctionUtility(rParameters["angular_velocity_constraints_settings"]["value"][i].GetString()));
            }
        }

        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyKinematicConstraintsProcess() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyKinematicConstraintsProcess algorithms.
    void Execute() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(!mInterval.IsInInterval(time)) return;

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++) {
                ModelPart::NodesContainerType::iterator it_node = it_begin + i;

                array_1d<double, 3>& vel = it_node->FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3>& ang_vel = it_node->FastGetSolutionStepValue(ANGULAR_VELOCITY);

                if(mVelocityIsConstrained[0]) {
                    it_node->Set(DEMFlags::FIXED_VEL_X, true);
                    it_node->pGetDof(VELOCITY_X)->FixDof();
                }
                if(mVelocityIsConstrained[1]) {
                    it_node->Set(DEMFlags::FIXED_VEL_Y, true);
                    it_node->pGetDof(VELOCITY_Y)->FixDof();
                }
                if(mVelocityIsConstrained[2]) {
                    it_node->Set(DEMFlags::FIXED_VEL_Z, true);
                    it_node->pGetDof(VELOCITY_Z)->FixDof();
                }
                if(mAngularVelocityIsConstrained[0]) {
                    it_node->Set(DEMFlags::FIXED_ANG_VEL_X, true);
                    it_node->pGetDof(ANGULAR_VELOCITY_X)->FixDof();
                }
                if(mAngularVelocityIsConstrained[1]) {
                    it_node->Set(DEMFlags::FIXED_ANG_VEL_Y, true);
                    it_node->pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
                }
                if(mAngularVelocityIsConstrained[2]) {
                    it_node->Set(DEMFlags::FIXED_ANG_VEL_Z, true);
                    it_node->pGetDof(ANGULAR_VELOCITY_Z)->FixDof();
                }

                for(int i=0; i<3; i++) {
                    if(mVelocityIsConstrained[i]) {
                        double velocity_value = 0.0;
                        if(mVelocityValueIsNumeric[i]) {
                            velocity_value = mVelocityValues[i];
                        }
                        else {
                            velocity_value = mVelocityFunctions[i].CallFunction(it_node->X(), it_node->Y(), it_node->Z(), time);
                        }
                        vel[i] = velocity_value;
                    }
                    if(mAngularVelocityIsConstrained[i]) {
                        double angular_velocity_value = 0.0;
                        if(mAngularVelocityValueIsNumeric[i]) {
                            angular_velocity_value = mAngularVelocityValues[i];
                        }
                        else {
                            angular_velocity_value = mAngularVelocityFunctions[i].CallFunction(it_node->X(), it_node->Y(), it_node->Z(), time);
                        }
                        ang_vel[i] = angular_velocity_value;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(!mInterval.IsInInterval(time)) return;

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++) {
                ModelPart::NodesContainerType::iterator it_node = it_begin + i;
                it_node->Set(DEMFlags::FIXED_VEL_X, false);
                it_node->Set(DEMFlags::FIXED_VEL_Y, false);
                it_node->Set(DEMFlags::FIXED_VEL_Z, false);
                it_node->Set(DEMFlags::FIXED_ANG_VEL_X, false);
                it_node->Set(DEMFlags::FIXED_ANG_VEL_Y, false);
                it_node->Set(DEMFlags::FIXED_ANG_VEL_Z, false);
                it_node->pGetDof(VELOCITY_X)->FreeDof();
                it_node->pGetDof(VELOCITY_Y)->FreeDof();
                it_node->pGetDof(VELOCITY_Z)->FreeDof();
                it_node->pGetDof(ANGULAR_VELOCITY_X)->FreeDof();
                it_node->pGetDof(ANGULAR_VELOCITY_Y)->FreeDof();
                it_node->pGetDof(ANGULAR_VELOCITY_Z)->FreeDof();
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyKinematicConstraintsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyKinematicConstraintsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    Parameters mParameters;
    IntervalUtility mInterval;
    array_1d<bool, 3> mVelocityIsConstrained;
    array_1d<bool, 3> mAngularVelocityIsConstrained;
    array_1d<bool, 3> mVelocityValueIsNumeric;
    array_1d<bool, 3> mAngularVelocityValueIsNumeric;
    array_1d<double, 3> mVelocityValues;
    array_1d<double, 3> mAngularVelocityValues;
    std::vector<PythonGenericFunctionUtility> mVelocityFunctions;
    std::vector<PythonGenericFunctionUtility> mAngularVelocityFunctions;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyKinematicConstraintsProcess& operator=(ApplyKinematicConstraintsProcess const& rOther);

    /// Copy constructor.
    //ApplyKinematicConstraintsProcess(ApplyKinematicConstraintsProcess const& rOther);

}; // Class ApplyKinematicConstraintsProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyKinematicConstraintsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyKinematicConstraintsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_KINEMATIC_CONSTRAINTS_PROCESS defined */
