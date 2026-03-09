//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_MANAGE_TIME_STEP_PROCESS_H_INCLUDED)
#define  KRATOS_MANAGE_TIME_STEP_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class ManageTimeStepProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ManageTimeStepProcess
    KRATOS_CLASS_POINTER_DEFINITION(ManageTimeStepProcess);

    ///@}
    ///@name Life Cycle
    ///@{


    ManageTimeStepProcess(ModelPart& rModelPart,
                          double& rMinDeltaTime, double& rMaxDeltaTime,
                          double& rReductionFactor, double& rIncreaseFactor,
                          double& rErrorTolerance, unsigned int& rMinIterations,
                          unsigned int& rMaxIterations, unsigned int& rNumberOfConstantSteps) : Process(Flags()),  mrModelPart(rModelPart), mMinDeltaTime(rMinDeltaTime), mMaxDeltaTime(rMaxDeltaTime), mReductionFactor(rReductionFactor), mIncreaseFactor(rIncreaseFactor), mErrorTolerance(rErrorTolerance), mMinIterations(rMinIterations), mMaxIterations(rMaxIterations), mNumberOfConstantSteps(rNumberOfConstantSteps)
    {
    }

    ManageTimeStepProcess( ModelPart& rModelPart,
                           Parameters CustomParameters )
      : mrModelPart(rModelPart)
    {
      KRATOS_TRY

      Parameters DefaultParameters( R"(
            {
                 "time_step": 1.0,
                 "start_time": 0.0,
                 "end_time": 1.0,
                 "adaptive_time_step":{
                     "minimum_time_step": 0.1,
                     "maximum_time_step": 1.0,
                     "reduction_factor": 2.0,
                     "increase_factor": 1.0,
                     "error_tolerance": 1e-4,
                     "minimum_iterations": 2,
                     "maximum_iterations": 10,
                     "number_constant_steps": 4
                 }
            }  )" );


      mAdaptiveTimeStep = false;
      if( CustomParameters.Has("adaptive_time_step") )
        mAdaptiveTimeStep = true;

      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);

      ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

      bool restarted = false;
      if( rCurrentProcessInfo.Has(IS_RESTARTED) ){
        if( rCurrentProcessInfo[IS_RESTARTED] == true ){
          restarted = true;
        }
      }

      if( !restarted ){
        rCurrentProcessInfo.SetValue(DELTA_TIME, CustomParameters["time_step"].GetDouble());
        rCurrentProcessInfo.SetValue(TIME, CustomParameters["start_time"].GetDouble());
      }

      mTime = rCurrentProcessInfo[TIME];
      mStep = rCurrentProcessInfo[STEP];
      mEndTime = CustomParameters["end_time"].GetDouble();

      if( mAdaptiveTimeStep )
        SetAdaptiveTimeParameters(CustomParameters["adaptive_time_step"]);

      KRATOS_CATCH(" ")

    }

    /// Destructor.
    virtual ~ManageTimeStepProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the ManageTimeStepProcess algorithms.
    void Execute() override
    {
        KRATOS_TRY;

//         const int nelements = mrModelPart.Elements().size();

//         if (nelements != 0)
//         {
//           ModelPart::ElementsContainerType::iterator it_begin =
//               mrModelPart.ElementsBegin();

// #pragma omp parallel for
//           for (int i = 0; i < nelements; i++)
//           {
//             ModelPart::ElementsContainerType::iterator it = it_begin + i;

//             if (this->MatchTransferFlags(*(it.base())))
//             {
//               this->AssignFlags(*(it.base()));
//             }
//           }
//         }



        KRATOS_CATCH("");
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
      ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

      if( mAdaptiveTimeStep )
        PredictTimeStep();

      double mTime = rCurrentProcessInfo[TIME] + rCurrentProcessInfo[DELTA_TIME];

      mrModelPart.ProcessInfo[STEP] = (++mStep);

      mrModelPart.CloneTimeStep(mTime);
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
    {
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ManageTimeStepProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ManageTimeStepProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    ManageTimeStepProcess(ManageTimeStepProcess const& rOther);

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

    bool mAdaptiveTimStep;

    double mTime;
    double mStep;
    double mEndTime;

    double mMinDeltaTime;
    double mMaxDeltaTime;

    double mReductionFactor;
    double mIncreaseFactor;

    double mErrorTolerance;

    unsigned int mMinIterations;
    unsigned int mMaxIterations;
    unsigned int mNumberOfConstantSteps;

    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{

    void SetAdaptiveTimeParameters(Parameters CustomParameters)
    {
      KRATOS_TRY

      Parameters DefaultParameters( R"(
      {
            "minimum_time_step": 0.1,
            "maximum_time_step": 1.0,
            "reduction_factor": 2.0,
            "increase_factor": 1.0,
            "error_tolerance": 1e-4,
            "minimum_iterations": 2,
            "maximum_iterations": 10,
            "number_constant_steps": 4
      }  )" );

      //validate against defaults -- this also ensures no type mismatch
      CustomParameters.ValidateAndAssignDefaults(DefaultParameters);


      mMinDeltaTime = CustomParameters["minimum_time_step"].GetDouble();
      mMaxDeltaTime = CustomParameters["maximum_time_step"].GetDouble();

      mReductionFactor = CustomParameters["reduction_factor"].GetDouble();
      mIncreaseFactor = CustomParameters["increase_factor"].GetDouble();

      mErrorTolerance = CustomParameters["error_tolerance"].GetDouble();

      mMinIterations = CustomParameters["minimum_iterations"].GetInt();
      mMaxIterations = CustomParameters["maximum_iterations"].GetInt();
      mNumberOfConstantSteps = CustomParameters["number_constant_steps"].GetInt();

      KRATOS_CATCH(" ")
    }
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    ManageTimeStepProcess& operator=(ManageTimeStepProcess const& rOther);

    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ManageTimeStepProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ManageTimeStepProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ManageTimeStepProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MANAGE_TIME_STEP_PROCESS_H_INCLUDED  defined
