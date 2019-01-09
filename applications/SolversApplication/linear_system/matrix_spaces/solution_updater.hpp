//
//   Project Name:        KratosSolversApplication    $
//   Created by:          $Author:        JMCarbonell $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:         January 2019 $
//   Revision:            $Revision:              0.0 $
//
//

#if !defined(KRATOS_SOLUTION_UPDATER_H_INCLUDED)
#define KRATOS_SOLUTION_UPDATER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup SolversApplication
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

/// Utility class to update the values of degree of freedom (Solution) variables after solving the system.
/** This class encapsulates the operation of updating nodal degrees of freedom after a system solution.
 *  In pseudo-code, the operation to be performed is
 *  for each dof: dof.variable += dx[dof.equation_id]
 *  This operation is a simple loop in shared memory, but requires additional infrastructure in distributed
 *  memory (MPI) in order to obtain out-of-process update data.
 *  SolutionUpdater takes care of both the operation and the eventual auxiliary infrastructure.
 */
template< class TSparseSpace >
class SolutionUpdater
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of SolutionUpdater
  KRATOS_CLASS_POINTER_DEFINITION(SolutionUpdater);

  using DofType = Dof<typename TSparseSpace::DataType>;
  using DofsArrayType = PointerVectorSet<DofType,
                                         SetIdentityFunction<DofType>,
                                         std::less<typename SetIdentityFunction<DofType>::result_type>,
                                         std::equal_to<typename SetIdentityFunction<DofType>::result_type>,
                                         DofType*>;

  using SystemVectorType = typename TSparseSpace::VectorType;
  using DataType = typename TSparseSpace::DataType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  SolutionUpdater() {}

  /// Deleted copy constructor
  SolutionUpdater(SolutionUpdater const& rOther) = delete;

  /// Destructor.
  virtual ~SolutionUpdater() {}

  /// Deleted assignment operator
  SolutionUpdater& operator=(SolutionUpdater const& rOther) = delete;

  ///@}
  ///@name Operations
  ///@{

  /// Create a new instance of this class.
  /** This function is used by the SparseSpace class to create new
   *  SolutionUpdater instances of the appropriate type.
   *  @return a std::unique_pointer to the new instance.
   *  @see UblasSpace::CreateSolutionUpdater(), TrilinosSpace::CreateSolutionUpdater().
   */
  typename SolutionUpdater::UniquePointer Create() const
  {
    return Kratos::make_unique<SolutionUpdater>();
  }

  /// Initialize the SolutionUpdater in preparation for a subsequent UpdateSolution call.
  /** Note that the base SolutionUpdater does not have internal data, so this does nothing.
   *  @param[in] rDofSet The list of degrees of freedom.
   *  @param[in] rDx The update vector.
   */
  virtual void Initialize(const DofsArrayType& rDofSet, const SystemVectorType& rDx) {}

  /// Free internal storage to reset the instance and/or optimize memory consumption.
  /** Note that the base SolutionUpdater does not have internal data, so this does nothing.
   */
  virtual void Clear() {}

  /// Calculate new values for the problem's degrees of freedom using the update vector rDx.
  /** For each Dof in rDofSet, this function calculates the updated value for the corresponding
   *  variable as value += rDx[dof.EquationId()].
   *  @param[in/out] rDofSet The list of degrees of freedom.
   *  @param[in] rDx The update vector.
   *  This method will check if Initialize() was called before and call it if necessary.
   */
  virtual void AddSolution(DofsArrayType& rDofSet, const SystemVectorType& rDx)
  {
    const int num_dof = static_cast<int>(rDofSet.size());

#pragma omp parallel for
    for(int i = 0;  i < num_dof; ++i) {
      auto it_dof = rDofSet.begin() + i;

      if (it_dof->IsFree())
        it_dof->GetSolutionStepValue() += TSparseSpace::GetValue(rDx,it_dof->EquationId());
    }
  }

  virtual void SetSolution(DofsArrayType& rDofSet, const SystemVectorType& rDx)
  {
    const int num_dof = static_cast<int>(rDofSet.size());

#pragma omp parallel for
    for(int i = 0;  i < num_dof; ++i) {
      auto it_dof = rDofSet.begin() + i;

      if (it_dof->IsFree())
        it_dof->GetSolutionStepValue() = TSparseSpace::GetValue(rDx,it_dof->EquationId());
    }
  }


  virtual void AddSolution(std::vector<double*>& ValueArray, const std::vector<std::size_t>& IndexArray, const SystemVectorType& rDx)
  {
    const int index_size = static_cast<int>(IndexArray.size());

    DataType* values = new DataType(index_size);
    TSparseSpace::GatherValues(rDx, IndexArray, values);

#pragma omp parallel for
    for(int i = 0;  i < index_size; ++i)
    {
      *(ValueArray[i]) += values[i];
    }

    delete [] values;
  }

  virtual void SetSolution(std::vector<double*>& ValueArray, const std::vector<std::size_t>& IndexArray, const SystemVectorType& rDx)
  {
    const int index_size = static_cast<int>(IndexArray.size());

    DataType* values = new DataType(index_size);
    TSparseSpace::GatherValues(rDx, IndexArray, values);

#pragma omp parallel for
    for(int i = 0;  i < index_size; ++i)
    {
      *(ValueArray[i]) = values[i];
    }

    delete [] values;
  }

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const
  {
    std::stringstream buffer;
    buffer << "SolutionUpdater" ;
    return buffer.str();
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << this->Info() << std::endl;
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const
  {
    rOStream << this->Info() << std::endl;
  }

  ///@}
 protected:
  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Member Variables
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
  ///@name Serialization
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  ///@}

}; // Class SolutionUpdater


///@}
///@name Input and output
///@{

/// input stream function
template< class TSparseSpace >
inline std::istream& operator >> (std::istream& rIStream, SolutionUpdater<TSparseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TSparseSpace >
inline std::ostream& operator << (std::ostream& rOStream,const SolutionUpdater<TSparseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SOLUTION_UPDATER_H_INCLUDED  defined
