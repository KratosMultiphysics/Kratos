// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef RESPONSE_FUNCTION_H
#define RESPONSE_FUNCTION_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "shape_optimization_application.h"

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

class ResponseFunction
{
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of ResponseFunction
  KRATOS_CLASS_POINTER_DEFINITION(ResponseFunction);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  ResponseFunction()
  {
  }

  /// Destructor.
  virtual ~ResponseFunction()
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
  virtual void CalculateValue() = 0;

  // --------------------------------------------------------------------------
  virtual void CalculateGradient() = 0;

  // --------------------------------------------------------------------------
  virtual double GetValue() = 0;

  // --------------------------------------------------------------------------
  virtual double GetInitialValue() = 0;

  // --------------------------------------------------------------------------
  virtual boost::python::dict GetGradient() = 0;

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
    return "ResponseFunction";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "ResponseFunction";
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
  //      ResponseFunction& operator=(ResponseFunction const& rOther);

  /// Copy constructor.
  //      ResponseFunction(ResponseFunction const& rOther);

  ///@}

}; // Class ResponseFunction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RESPONSE_FUNCTION_H
