//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_CONVERGENCE_CRITERION_H_INCLUDED )
#define  KRATOS_CONVERGENCE_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_local_flags.hpp"
#include "includes/model_part.h"
//#include "includes/dof.h"

namespace Kratos
{
///@addtogroup SolidMechanicsApplication
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


/** @brief Convergence Criterion base class
 *  @details This is the base class for the convergence criteria
 */
template<class TSparseSpace, class TDenseSpace>
class ConvergenceCriterion : public Flags
{
 public:

  ///@name Type Definitions
  ///@{
  typedef CriterionLocalFlags                       LocalFlagType;
  typedef typename TSparseSpace::DataType                DataType;
  typedef ModelPart::DofsArrayType                  DofsArrayType;
  typedef typename TSparseSpace::MatrixType      SystemMatrixType;
  typedef typename TSparseSpace::VectorType      SystemVectorType;
  typedef typename TDenseSpace::MatrixType  LocalSystemMatrixType;
  typedef typename TDenseSpace::VectorType  LocalSystemVectorType;

  /// Pointer definition of ConvergenceCriterion
  KRATOS_CLASS_POINTER_DEFINITION(ConvergenceCriterion);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Constructor
  ConvergenceCriterion()
  {
    SetEchoLevel(1);
  }

  /// Copy contructor.
  ConvergenceCriterion( ConvergenceCriterion const& rOther)
      :mEchoLevel(rOther.mEchoLevel)
  {
  }

  /// Destructor.
  ~ConvergenceCriterion() override
  {
  }

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  //*********************************************************************************

  /**level of echo for the convergence criterion
     0 -> mute... no echo at all
     1 -> print basic informations
     2 -> print extra informations
  */
  virtual void SetEchoLevel(int Level)
  {
    mEchoLevel = Level;
  }

  int GetEchoLevel()
  {
    return mEchoLevel;
  }


  /*Criterias that need to be called before getting the solution */
  virtual bool PreCriteria(ModelPart& rModelPart,
                           DofsArrayType& rDofSet,
                           const SystemMatrixType& rA,
                           const SystemVectorType& rDx,
                           const SystemVectorType& rb)
  {
    return true;
  }

  /*Criterias that need to be called after getting the solution */
  virtual bool PostCriteria(ModelPart& rModelPart,
                            DofsArrayType& rDofSet,
                            const SystemMatrixType& rA,
                            const SystemVectorType& rDx,
                            const SystemVectorType& rb)

  {
    return true;
  }

  virtual void InitializeSolutionStep(ModelPart& rModelPart,
                                      DofsArrayType& rDofSet,
                                      const SystemMatrixType& rA,
                                      const SystemVectorType& rDx,
                                      const SystemVectorType& rb)
  {
  }

  virtual void FinalizeSolutionStep(ModelPart& rModelPart,
                                    DofsArrayType& rDofSet,
                                    const SystemMatrixType& rA,
                                    const SystemVectorType& rDx,
                                    const SystemVectorType& rb)
  {
  }

  /**
   * This function is designed to be called once to perform all the checks needed
   * on the input provided. Checks can be "expensive" as the function is designed
   * to catch user's errors.
   * @param rModelPart
   * @return 0 all ok
   */
  virtual int Check(ModelPart& rModelPart)
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH("")
  }

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
  ///@name Friends
  ///@{

  ///@}

 protected:
  ///@name Protected static Member Variables
  ///@{

  ///@}
  ///@name Protected member Variables
  ///@{

  int  mEchoLevel;

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

  ///@}

};  // Class ConvergenceCriterion

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CONVERGENCE_CRITERION_H_INCLUDED defined
