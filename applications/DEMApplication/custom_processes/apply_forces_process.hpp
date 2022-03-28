//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Joaquín Irazábal González
//

#if !defined(KRATOS_APPLY_FORCES_PROCESS )
#define  KRATOS_APPLY_FORCES_PROCESS


// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/table.h"
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


  /** This process sets the EXTERNAL_APPLIED_FORCE variables.
   */
  class KRATOS_API(DEM_APPLICATION) ApplyForcesProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of ApplyForcesProcess
      KRATOS_CLASS_POINTER_DEFINITION(ApplyForcesProcess);

      /// Defining a table with double argument and result type as table type.
      typedef Table<double,double> TableType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor
      ApplyForcesProcess(ModelPart& rModelPart, Parameters rParameters);

      /// Destructor.
      ~ApplyForcesProcess() override;


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Execute() override;

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

      ModelPart& mrModelPart;
      Parameters mParameters;
      IntervalUtility mInterval;
      array_1d<bool, 3> mForceValueIsNumeric;
      array_1d<double, 3> mForceValues;
      std::vector<GenericFunctionUtility> mForceFunctions;
      array_1d<int, 3> mForceTableId;
      std::vector<TableType::Pointer> mpForceTable;

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
      ApplyForcesProcess& operator=(ApplyForcesProcess const& rOther);

      /// Copy constructor.
      ApplyForcesProcess(ApplyForcesProcess const& rOther);

      ///@}

    }; // Class ApplyForcesProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const ApplyForcesProcess& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_APPLY_FORCES_PROCESS  defined
