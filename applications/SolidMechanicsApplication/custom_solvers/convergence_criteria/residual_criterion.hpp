//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_CRITERION_H_INCLUDED )
#define  KRATOS_RESIDUAL_CRITERION_H_INCLUDED

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
 * @class ResidualCriterion
 * @brief This convergence criteria checks the residual
 */
template<class TSparseSpace, class TDenseSpace>
class ResidualCriterion : public  ConvergenceCriterion< TSparseSpace, TDenseSpace >
{
 public:

  ///@name Type Definitions
  ///@{

  typedef ConvergenceCriterion< TSparseSpace, TDenseSpace > BaseType;
  typedef typename BaseType::LocalFlagType             LocalFlagType;
  typedef typename BaseType::DataType                       DataType;
  typedef typename BaseType::DofsArrayType             DofsArrayType;
  typedef typename BaseType::SystemMatrixType       SystemMatrixType;
  typedef typename BaseType::SystemVectorType       SystemVectorType;

  /// Pointer definition of ResidualCriterion
  KRATOS_CLASS_POINTER_DEFINITION(ResidualCriterion);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Constructor.
  ResidualCriterion(DataType RatioTolerance,
                    DataType AlwaysConvergedNorm)
      : BaseType()
  {
    mRatioTolerance       = RatioTolerance;
    mAlwaysConvergedNorm  = AlwaysConvergedNorm;
  }

  /// Copy constructor.
  ResidualCriterion( ResidualCriterion const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mInitialResidualNorm(rOther.mInitialResidualNorm)
      ,mCurrentResidualNorm(rOther.mCurrentResidualNorm)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
  {
  }

  /// Destructor.
  ~ResidualCriterion() override {}


  ///@}
  ///@name Operators
  ///@{

  //Criterion that needs to be called after getting the solution
  bool PostCriteria(ModelPart& rModelPart,
                    DofsArrayType& rDofSet,
                    const SystemMatrixType& rA,
                    const SystemVectorType& rDx,
                    const SystemVectorType& rb
                    ) override
  {
    if (TSparseSpace::Size(rb) != 0) //if we are solving for something
    {

      if( this->IsNot(LocalFlagType::INITIALIZED) )
      {
        mInitialResidualNorm = TSparseSpace::TwoNorm(rb);
        this->Set(LocalFlagType::INITIALIZED,true);
      }

      DataType ratio;
      mCurrentResidualNorm = TSparseSpace::TwoNorm(rb);

      if(mInitialResidualNorm == 0.00)
      {
        ratio = 0.00;
      }

      else
      {
        ratio = mCurrentResidualNorm/mInitialResidualNorm;
      }

      double b_size = TSparseSpace::Size(rb);
      DataType absolute_norm = (mCurrentResidualNorm/b_size);

      if (rModelPart.GetCommunicator().MyPID() == 0)
      {
        if (this->GetEchoLevel() >= 1)
        {
          std::cout << "RESIDUAL CRITERION :: Ratio = "<< ratio  << ";  Norm = " << absolute_norm << std::endl;
        }
      }

      rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
      rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

      if (ratio <= mRatioTolerance || absolute_norm < mAlwaysConvergedNorm)
      {
        if (rModelPart.GetCommunicator().MyPID() == 0)
        {
          if (this->GetEchoLevel() >= 1)
          {
            std::cout << "Convergence is achieved" << std::endl;
          }
        }
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return true;
    }
  }


  void InitializeSolutionStep(ModelPart& rModelPart,
                              DofsArrayType& rDofSet,
                              const SystemMatrixType& rA,
                              const SystemVectorType& rDx,
                              const SystemVectorType& rb) override
  {
    this->Set(LocalFlagType::INITIALIZED,false);
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

  DataType mRatioTolerance;

  DataType mInitialResidualNorm;

  DataType mCurrentResidualNorm;

  DataType mAlwaysConvergedNorm;

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

}; // Class ResidualCriterion

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RESIDUAL_CRITERION_H_INCLUDED defined
