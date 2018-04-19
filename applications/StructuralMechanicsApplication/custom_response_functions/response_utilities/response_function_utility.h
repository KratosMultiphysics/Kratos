// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//	             license: structural_mechanics_application/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//

#ifndef RESPONSE_FUNCTION_UTILITY_H
#define RESPONSE_FUNCTION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"

// ==============================================================================

namespace Kratos
{

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

/// Short class definition.
/** Detail class definition.

*/

class ResponseFunctionUtility
{
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of ResponseFunctionUtility
  KRATOS_CLASS_POINTER_DEFINITION(ResponseFunctionUtility);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  ResponseFunctionUtility()
  {
  }

  /// Destructor.
  virtual ~ResponseFunctionUtility()
  {
  }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  // ==============================================================================
  virtual void Initialize() = 0;

  // --------------------------------------------------------------------------
  virtual double CalculateValue() = 0;

  // --------------------------------------------------------------------------
  virtual void CalculateGradient() = 0;

  // --------------------------------------------------------------------------
  virtual void ConsiderDiscretization() = 0;

  // ==============================================================================

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
  virtual std::string Info() const
  {
    return "ResponseFunctionUtility";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "ResponseFunctionUtility";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream &rOStream) const
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
  //      ResponseFunctionUtility& operator=(ResponseFunctionUtility const& rOther);

  /// Copy constructor.
  //      ResponseFunctionUtility(ResponseFunctionUtility const& rOther);

  ///@}

}; // Class ResponseFunctionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RESPONSE_FUNCTION_H
