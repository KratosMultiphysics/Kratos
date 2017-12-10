//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_REGISTER_LINEAR_SOLVERS )
#define  KRATOS_REGISTER_LINEAR_SOLVERS



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

namespace Kratos
{
  ///@addtogroup ExternalSolversApplication
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

  /// registers the linear solvers to kratos
  /** registers the linear solvers to kratos
  */
  class RegisterLinearSolvers
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of RegisterLinearSolvers
      KRATOS_CLASS_POINTER_DEFINITION(RegisterLinearSolvers);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      RegisterLinearSolvers();

      /// Destructor.
      virtual ~RegisterLinearSolvers(){};


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{


      ///@}
      ///@name Access
      ///@{


      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{

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
      RegisterLinearSolvers& operator=(RegisterLinearSolvers const& rOther) = delete;

      /// Copy constructor.
      RegisterLinearSolvers(RegisterLinearSolvers const& rOther) = delete;


      ///@}

    }; // Class RegisterLinearSolvers

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_REGISTER_LINEAR_SOLVERS