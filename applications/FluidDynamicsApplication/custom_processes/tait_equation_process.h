//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#if !defined(KRATOS_TAIT_EQUATION_PROCESS_H_INCLUDED)
#define  KRATOS_TAIT_EQUATION_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
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

 
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TaitEquationProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TaitEquationProcess
      KRATOS_CLASS_POINTER_DEFINITION(TaitEquationProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor
      TaitEquationProcess(ModelPart& rModelPart, Parameters& rParameters);

      /// Destructor.
      ~TaitEquationProcess() override;


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Execute() override;

      void ExecuteInitialize() override;

      void ExecuteInitializeSolutionStep() override;

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

      void ValidateModelPart();

      void AssignTaitEquation();

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

      double mtheta, mk_0, mp_0, mrho_0;

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
      TaitEquationProcess& operator=(TaitEquationProcess const& rOther);

      /// Copy constructor.
      TaitEquationProcess(TaitEquationProcess const& rOther);


      ///@}

    }; // Class TaitEquationProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const TaitEquationProcess& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUSSINESQ_FORCE_PROCESS_H_INCLUDED  defined
