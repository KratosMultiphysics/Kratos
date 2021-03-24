//
//   Author:   Joaquín Irazábal González
//

#if !defined(KRATOS_APPLY_FORCES_AND_MOMENTS_PROCESS )
#define  KRATOS_APPLY_FORCES_AND_MOMENTS_PROCESS

#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/interval_utility.h"
#include "utilities/python_function_callback_utility.h"
#include "DEM_application_variables.h"

namespace Kratos
{

class ApplyForcesAndMomentsProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyForcesAndMomentsProcess);


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyForcesAndMomentsProcess(ModelPart& model_part,
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
                "force_settings" : {
                    "value"                : [10.0, "3*t", "x+y"]
                },
                "moment_settings" : {
                    "value"                : [10.0, "3*t", "x+y"]
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

        for(int i=0; i<3; i++) {
            if(rParameters["force_settings"]["value"][i].IsNumber()) {
                mForceValueIsNumeric[i] = true;
                mForceValues[i] = rParameters["force_settings"]["value"][i].GetDouble();
                mForceFunctions.push_back(PythonGenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                mForceValueIsNumeric[i] = false;
                mForceFunctions.push_back(PythonGenericFunctionUtility(rParameters["force_settings"]["value"][i].GetString()));
            }
            if(rParameters["moment_settings"]["value"][i].IsNumber()) {
                mMomentValueIsNumeric[i] = true;
                mMomentValues[i] = rParameters["moment_settings"]["value"][i].GetDouble();
                mMomentFunctions.push_back(PythonGenericFunctionUtility("0.0")); // because I can't construct an array_1d of these
            }
            else {
                mMomentValueIsNumeric[i] = false;
                mMomentFunctions.push_back(PythonGenericFunctionUtility(rParameters["moment_settings"]["value"][i].GetString()));
            }
        }

        mParameters = rParameters;

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyForcesAndMomentsProcess() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyForcesAndMomentsProcess algorithms.
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

                array_1d<double, 3>& force = it_node->FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
                array_1d<double, 3>& moment = it_node->FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);

                for(int i=0; i<3; i++) {
                    double force_value = 0.0;
                    if(mForceValueIsNumeric[i]) {
                        force_value = mForceValues[i];
                    }
                    else {
                        force_value = mForceFunctions[i].CallFunction(it_node->X(), it_node->Y(), it_node->Z(), time);
                    }
                    force[i] = force_value;

                    double moment_value = 0.0;
                    if(mMomentValueIsNumeric[i]) {
                        moment_value = mMomentValues[i];
                    }
                    else {
                            moment_value = mMomentFunctions[i].CallFunction(it_node->X(), it_node->Y(), it_node->Z(), time);
                    }
                    moment[i] = moment_value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;


        const double time = mrModelPart.GetProcessInfo()[TIME];

        if(mInterval.IsInInterval(time)) return;

        const int nnodes = static_cast<int>(mrModelPart.Nodes().size());

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++) {
                ModelPart::NodesContainerType::iterator it_node = it_begin + i;
                it_node->FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE) = ZeroVector(3);
                it_node->FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT) = ZeroVector(3);
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyForcesAndMomentsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyForcesAndMomentsProcess";
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
    array_1d<bool, 3> mForceValueIsNumeric;
    array_1d<bool, 3> mMomentValueIsNumeric;
    array_1d<double, 3> mForceValues;
    array_1d<double, 3> mMomentValues;
    std::vector<PythonGenericFunctionUtility> mForceFunctions;
    std::vector<PythonGenericFunctionUtility> mMomentFunctions;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyForcesAndMomentsProcess& operator=(ApplyForcesAndMomentsProcess const& rOther);

    /// Copy constructor.
    //ApplyForcesAndMomentsProcess(ApplyForcesAndMomentsProcess const& rOther);

}; // Class ApplyForcesAndMomentsProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyForcesAndMomentsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyForcesAndMomentsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_FORCES_AND_MOMENTS_PROCESS defined */
