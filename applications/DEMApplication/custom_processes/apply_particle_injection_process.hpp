
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt



#if !defined(KRATOS_APPLY_PARTICLE_INJECTION_PROCESS )
#define  KRATOS_APPLY_PARTICLE_INJECTION_PROCESS

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/interval_utility.h"
#include "utilities/python_function_callback_utility.h"
#include "DEM_application_variables.h"

namespace Kratos
{
  ///@addtogroup DEMApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{


  class KRATOS_API(DEM_APPLICATION) ApplyParticleInjectionProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      KRATOS_CLASS_POINTER_DEFINITION(ApplyParticleInjectionProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      ApplyParticleInjectionProcess(ModelPart& rInletPart, ModelPart& rSpheresPart, Parameters rParameters);

      /// Destructor.
      ~ApplyParticleInjectionProcess() {};


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Execute() override;

      void ExecuteInitialize() override;

      void ExecuteInitializeSolutionStep() override;

      void ExecuteFinalizeSolutionStep() override;

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
      std::string Info() const override;

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override;


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

      ModelPart& mrInletPart;
      ModelPart& mrSpheresPart;
      Parameters mParameters;
      IntervalUtility mInterval;  // this could be used instead of inlet_start_time and stop_time
      bool mStrategyForContinuum;

      // array_1d<bool, 3> mVelocityIsConstrained;
      // array_1d<bool, 3> mVelocityValueIsNumeric;
      // array_1d<double, 3> mVelocityValues;
      // std::vector<GenericFunctionUtility> mVelocityFunctions;


      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{


      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      ApplyParticleInjectionProcess& operator=(ApplyParticleInjectionProcess const& rOther);

      /// Copy constructor.
      ApplyParticleInjectionProcess(ApplyParticleInjectionProcess const& rOther);


      ///@}

    };

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const ApplyParticleInjectionProcess& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif
