// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef RESPONSE_H
#define RESPONSE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "optimization_application.h"

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

class Response
{
public:
  ///@name Type Definitions
  ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;

  /// Pointer definition of Response
  KRATOS_CLASS_POINTER_DEFINITION(Response);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  Response(std::string ResponseName, std::string ResponseType, Parameters ResponseSettings): mResponseName(ResponseName), mResponseType(ResponseType), mResponseSettings(ResponseSettings){}

  /// Destructor.
  virtual ~Response(){}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  // --------------------------------------------------------------------------
  virtual void Initialize(){};

  // --------------------------------------------------------------------------
  virtual double CalculateValue(){};

  // --------------------------------------------------------------------------
  virtual void CalculateGradient() {};

  // --------------------------------------------------------------------------

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
    return "Response base class";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "Response base class";
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
    std::string mResponseName;
    std::string mResponseType;
    Parameters mResponseSettings;

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
  //      RESPONSE& operator=(RESPONSE const& rOther);

  /// Copy constructor.
  //      RESPONSE(RESPONSE const& rOther);

  ///@}

}; // Class RESPONSE

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RESPONSE_H
