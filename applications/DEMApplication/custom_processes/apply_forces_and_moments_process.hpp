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

#if !defined(KRATOS_APPLY_FORCES_AND_MOMENTS_PROCESS )
#define  KRATOS_APPLY_FORCES_AND_MOMENTS_PROCESS



// System includes
#include <string>
#include <iostream>


// External includes


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

  /// Auxiliary process to apply forces and moments to particles.
  /** This process sets the EXTERNAL_APPLIED_FORCE and EXTERNAL_APPLIED_MOMENT variables
      over particles.
   */
  class KRATOS_API(DEM_APPLICATION) ApplyForcesAndMomentsProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of ApplyForcesAndMomentsProcess
      KRATOS_CLASS_POINTER_DEFINITION(ApplyForcesAndMomentsProcess);

      /// Defining a table with double argument and result type as table type.
      typedef Table<double,double> TableType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor
      ApplyForcesAndMomentsProcess(ModelPart& rModelPart, Parameters rParameters);

      /// Destructor.
      ~ApplyForcesAndMomentsProcess() override;


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
      array_1d<bool, 3> mMomentValueIsNumeric;
      array_1d<double, 3> mForceValues;
      array_1d<double, 3> mMomentValues;
      std::vector<GenericFunctionUtility> mForceFunctions;
      std::vector<GenericFunctionUtility> mMomentFunctions;
      array_1d<int, 3> mForceTableId;
      array_1d<int, 3> mMomentTableId;
      std::vector<TableType::Pointer> mpForceTable;
      std::vector<TableType::Pointer> mpMomentTable;

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
      ApplyForcesAndMomentsProcess& operator=(ApplyForcesAndMomentsProcess const& rOther);

      /// Copy constructor.
      ApplyForcesAndMomentsProcess(ApplyForcesAndMomentsProcess const& rOther);


      ///@}

    }; // Class ApplyForcesAndMomentsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const ApplyForcesAndMomentsProcess& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_APPLY_FORCES_AND_MOMENTS_PROCESS  defined
