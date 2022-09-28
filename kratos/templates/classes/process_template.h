//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED )
#define  KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED


// System includes
@{KRATOS_SYSTEM_INCLUDES}
#include <string>
#include <iostream>
#include <algorithm>

// External includes
@{KRATOS_EXTERNA_INCLUDES}
#include "includes/kratos_flags.h"

// Project includes
@{KRATOS_PROJECT_INCLUDES}
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "elements/distance_calculation_element_simplex.h"

namespace Kratos {

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
@{KRATOS_CLASS_TEMPLATE}
class @{KRATOS_NAME_CAMEL} : public @{KRATOS_CLASS_BASE} {
public:

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of @{KRATOS_NAME_CAMEL}
  KRATOS_CLASS_POINTER_DEFINITION(@{KRATOS_NAME_CAMEL});

  ///@}
  ///@name Life Cycle
  ///@{

  @{KRATOS_NAME_CAMEL}() {
  }

  /// Destructor.
  virtual ~@{KRATOS_NAME_CAMEL}() {
  }

  ///@}
  ///@name Operators
  ///@{

  void operator()() {
    Execute();
  }

  ///@}
  ///@name Operations
  ///@{

  virtual void Execute() {
  }

  virtual void Clear() {
  }

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
  virtual std::string Info() const {
      return "@{KRATOS_NAME_CAMEL}";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "@{KRATOS_NAME_CAMEL}";
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const {
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
  @{KRATOS_NAME_CAMEL}& operator=(@{KRATOS_NAME_CAMEL} const& rOther);

  /// Copy constructor.
  //@{KRATOS_NAME_CAMEL}(@{KRATOS_NAME_CAMEL} const& rOther);

  ///@}

}; // Class @{KRATOS_NAME_CAMEL}

//avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags @{KRATOS_NAME_CAMEL}<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::PERFORM_STEP1(Kratos::Flags::Create(0));

template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags @{KRATOS_NAME_CAMEL}<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (std::istream& rIStream,
                                  @{KRATOS_NAME_CAMEL}<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis);

/// output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const @{KRATOS_NAME_CAMEL}<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED  defined
