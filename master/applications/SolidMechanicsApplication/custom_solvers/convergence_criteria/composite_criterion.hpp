//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPOSITE_CRITERION_H_INCLUDED )
#define  KRATOS_COMPOSITE_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/convergence_criteria/convergence_criterion.hpp"

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

/**
 * @class CompositeCriterion
 * @brief This convergence criteria checks simultaneously different convergence criteria
 */
template<class TSparseSpace,class TDenseSpace>
class CompositeCriterion : public ConvergenceCriterion< TSparseSpace, TDenseSpace >
{
 public:
  ///@name Type Definitions
  ///@{

  typedef ConvergenceCriterion<TSparseSpace,TDenseSpace>      BaseType;
  typedef typename BaseType::LocalFlagType               LocalFlagType;
  typedef typename BaseType::DataType                         DataType;
  typedef typename BaseType::DofsArrayType               DofsArrayType;
  typedef typename BaseType::SystemMatrixType         SystemMatrixType;
  typedef typename BaseType::SystemVectorType         SystemVectorType;

  typedef typename BaseType::Pointer                   BasePointerType;
  typedef std::vector<BasePointerType>  ConvergenceCriterionVectorType;

  /// Pointer definition of CompositeCriterion
  KRATOS_CLASS_POINTER_DEFINITION(CompositeCriterion);


  ///@}
  ///@name Life Cycle
  ///@{

  /// Constructor.
  CompositeCriterion(BasePointerType pFirstCriterion, BasePointerType pSecondCriterion) : BaseType()
  {
    mConvergenceCriteria.push_back(pFirstCriterion);
    mConvergenceCriteria.push_back(pSecondCriterion);
  }

  /// Constructor.
  CompositeCriterion(ConvergenceCriterionVectorType& rConvergenceCriteria) : BaseType(), mConvergenceCriteria(rConvergenceCriteria)
  {
  }

  /// Copy constructor
  CompositeCriterion(CompositeCriterion const& rOther)
      :BaseType(rOther)
      ,mConvergenceCriteria(rOther.mConvergenceCriteria)
  {
  }

  /// Destructor.
  ~CompositeCriterion () override {}

  ///@}
  ///@name Operators
  ///@{

  /**
   * @brief It sets the level of echo for the solving strategy
   * @param Level The level to set
   * @details The different levels of echo are:
   * - 0: Mute... no echo at all
   * - 1: Printing time and basic informations
   * - 2: Printing extra informations
   */
  void SetEchoLevel(int Level) override
  {
    BaseType::SetEchoLevel(Level);
    for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
        it!=mConvergenceCriteria.end(); ++it)
      (*it)->SetEchoLevel(Level);
  }

  /**
   * @return true if convergence is achieved, false otherwise
   */
  bool PreCriteria(ModelPart& rModelPart,
                   DofsArrayType& rDofSet,
                   const SystemMatrixType& rA,
                   const SystemVectorType& rDx,
                   const SystemVectorType& rb) override
  {
    if(this->Is(LocalFlagType::AND)){

      for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
          it!=mConvergenceCriteria.end(); ++it)
        if( !((*it)->PreCriteria(rModelPart,rDofSet,rA,rDx,rb)) )
          return false;

      return true;
    }
    else if (this->Is(LocalFlagType::OR)){

      for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
          it!=mConvergenceCriteria.end(); ++it)
        if( ((*it)->PreCriteria(rModelPart,rDofSet,rA,rDx,rb)) )
          return true;

      return false;
    }

    return false;
  }

  /**
   * @return true if convergence is achieved, false otherwise
   */
  bool PostCriteria(ModelPart& rModelPart,
                    DofsArrayType& rDofSet,
                    const SystemMatrixType& rA,
                    const SystemVectorType& rDx,
                    const SystemVectorType& rb) override
  {
    if(this->Is(LocalFlagType::AND)){

      for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
          it!=mConvergenceCriteria.end(); ++it)
        if( !((*it)->PostCriteria(rModelPart,rDofSet,rA,rDx,rb)) )
          return false;

      return true;
    }
    else if (this->Is(LocalFlagType::OR)){

      for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
          it!=mConvergenceCriteria.end(); ++it)
        if( ((*it)->PostCriteria(rModelPart,rDofSet,rA,rDx,rb)) )
          return true;

      return false;
    }

    return false;
  }

  void InitializeSolutionStep(ModelPart& rModelPart,
                              DofsArrayType& rDofSet,
                              const SystemMatrixType& rA,
                              const SystemVectorType& rDx,
                              const SystemVectorType& rb) override
  {
    for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
        it!=mConvergenceCriteria.end(); ++it)
      (*it)->InitializeSolutionStep(rModelPart,rDofSet,rA,rDx,rb);
  }

  void FinalizeSolutionStep(ModelPart& rModelPart,
                            DofsArrayType& rDofSet,
                            const SystemMatrixType& rA,
                            const SystemVectorType& rDx,
                            const SystemVectorType& rb) override
  {
    for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
        it!=mConvergenceCriteria.end(); ++it)
      (*it)->FinalizeSolutionStep(rModelPart,rDofSet,rA,rDx,rb);
  }


  /**
   * @brief This function is designed to be called once to perform all the checks needed on the input provided.
   * @details Checks can be "expensive" as the function is designed
   * to catch user's errors.
   * @param rModelPart ModelPart containing the problem.
   * @return 0 all ok
   */
  int Check(ModelPart& rModelPart) override
  {
    KRATOS_TRY

    if( !(this->Is(LocalFlagType::AND) || this->Is(LocalFlagType::OR)) )
      KRATOS_ERROR << "Flag AND or OR is not set for the Composite Convergence Criterion "<< std::endl;

    for(typename ConvergenceCriterionVectorType::iterator it=mConvergenceCriteria.begin();
        it!=mConvergenceCriteria.end(); ++it)
      if( ((*it)->Check(rModelPart)) != 0 )
        return 1;

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

  ConvergenceCriterionVectorType mConvergenceCriteria;

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

};  // Class CompositeCriterion

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COMPOSITE_CRITERION_H_INCLUDED  defined
