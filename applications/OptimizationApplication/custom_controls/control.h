// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef CONTROL_H
#define CONTROL_H

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

class Control
{
public:
  ///@name Type Definitions
  ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;

  /// Pointer definition of Control
  KRATOS_CLASS_POINTER_DEFINITION(Control);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  Control(std::string ControlName, std::string ControlType, Parameters ControlSettings): mControlName(ControlName), mControlType(ControlType), mControlSettings(ControlSettings){}

  /// Destructor.
  virtual ~Control(){}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  // --------------------------------------------------------------------------
  virtual void Initialize(){};

  // --------------------------------------------------------------------------
  virtual void Update(){};

  // --------------------------------------------------------------------------
  virtual void MapControlUpdate(const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable){};

  // --------------------------------------------------------------------------
  virtual void MapControlUpdate(const Variable<double> &rOriginVariable, const Variable<double> &rDestinationVariable) {};

  // --------------------------------------------------------------------------
  virtual void MapFirstDerivative(const Variable<array_3d> &rDestinationVariable, const Variable<array_3d> &rOriginVariable){};

  // --------------------------------------------------------------------------
  virtual void MapFirstDerivative(const Variable<double> &rDestinationVariable, const Variable<double> &rOriginVariable) {};

  // --------------------------------------------------------------------------
  virtual void Finalize(){};

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
    return "Control base class";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "Control base class";
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
    std::string mControlName;
    std::string mControlType;
    Parameters mControlSettings;

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
  //      CONTROL& operator=(CONTROL const& rOther);

  /// Copy constructor.
  //      CONTROL(CONTROL const& rOther);

  ///@}

}; // Class CONTROL

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // CONTROL_H
