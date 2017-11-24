//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED )
#define  KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED


// System includes
// Please put system includes in @{KRATOS_NAME_LOWER}.h

// External includes
// Please put external includes in @{KRATOS_NAME_LOWER}.h

// Project includes
@{KRATOS_PROJECT_INCLUDES}
#include "custom_processes/@{KRATOS_NAME_LOWER}.h"

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
class PostprocessEigenvaluesProcess : public @{KRATOS_CLASS_BASE} 
{
public:
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of PostprocessEigenvaluesProcess
  KRATOS_CLASS_POINTER_DEFINITION(PostprocessEigenvaluesProcess);

  ///@}
  ///@name Life Cycle
  ///@{

  PostprocessEigenvaluesProcess() {
  }

  /// Destructor.
  virtual ~PostprocessEigenvaluesProcess() {
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
      return "PostprocessEigenvaluesProcess";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "PostprocessEigenvaluesProcess";
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
  PostprocessEigenvaluesProcess& operator=(PostprocessEigenvaluesProcess const& rOther);

  /// Copy constructor.
  //PostprocessEigenvaluesProcess(PostprocessEigenvaluesProcess const& rOther);

  ///@}

}; // Class PostprocessEigenvaluesProcess

//avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags PostprocessEigenvaluesProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::PERFORM_STEP1(Kratos::Flags::Create(0));

template< unsigned int TDim,class TSparseSpace, class TDenseSpace, class TLinearSolver >
const Kratos::Flags PostprocessEigenvaluesProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (std::istream& rIStream,
                                  PostprocessEigenvaluesProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis);

/// output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PostprocessEigenvaluesProcess<TDim,TSparseSpace,TDenseSpace,TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED  defined
