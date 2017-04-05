// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef VTK_DATA_WRITER_UTILITIES_H
#define VTK_DATA_WRITER_UTILITIES_H

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

class OptimizationDataWriter
{
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of OptimizationDataWriter
  KRATOS_CLASS_POINTER_DEFINITION(OptimizationDataWriter);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  OptimizationDataWriter()
  {
  }

  /// Destructor.
  virtual ~OptimizationDataWriter()
  {
  }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  // ==============================================================================
  void initialize();


  // --------------------------------------------------------------------------
  virtual void initializeDesignOutput() = 0;

  // --------------------------------------------------------------------------
  virtual void initializeOptimizationLog() = 0;

  // --------------------------------------------------------------------------
  virtual double writeCurrentDesign() = 0;

  // --------------------------------------------------------------------------
  virtual double logCurrentOptimizationStep() = 0;

  // --------------------------------------------------------------------------
  virtual double finalizeDesignOutput() = 0;

  // --------------------------------------------------------------------------
  virtual double fianlizeOptimizationLog() = 0;

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
    return "OptimizationDataWriter";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream &rOStream) const
  {
    rOStream << "OptimizationDataWriter";
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
  //      OptimizationDataWriter& operator=(OptimizationDataWriter const& rOther);

  /// Copy constructor.
  //      OptimizationDataWriter(OptimizationDataWriter const& rOther);

  ///@}

}; // Class OptimizationDataWriter

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // VTK_DATA_WRITER_UTILITIES_H
