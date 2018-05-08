//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_VARIABLE_CRITERION_H_INCLUDED )
#define  KRATOS_VARIABLE_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/convergence_criteria/convergence_criterion.hpp"

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


/** @class VariableCriterion
 * @brief This convergence criteria checks the variable
 */

template<class TSparseSpace, class TDenseSpace>
class VariableCriterion : public ConvergenceCriterion<TSparseSpace,TDenseSpace>
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

  /// Pointer definition of VariableCriterion
  KRATOS_CLASS_POINTER_DEFINITION( VariableCriterion );

  ///@}
  ///@name Life Cycle
  ///@{

  /** Constructor.
   */
  VariableCriterion(DataType RatioTolerance,
                    DataType AlwaysConvergedNorm)
      : BaseType()
  {
    mRatioTolerance = RatioTolerance;
    mAlwaysConvergedNorm = AlwaysConvergedNorm;
  }

  /// Copy constructor.
    
  VariableCriterion(VariableCriterion const& rOther)
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAlwaysConvergedNorm(rOther.mAlwaysConvergedNorm)
  {
  }

  /// Destructor.
  ~VariableCriterion() override {}


  ///@}
  ///@name Operators
  ///@{

  //Criterion that needs to be called after getting the solution
  bool PostCriteria(ModelPart& rModelPart,
                    DofsArrayType& rDofSet,
                    const SystemMatrixType& rA,
                    const SystemVectorType& rDx,
                    const SystemVectorType& rb) override
  {
    if (TSparseSpace::Size(rDx) != 0) //if we are solving for something
    {
      DataType ratio = 0.00;
      DataType CorrectionNorm = DataType();
      DataType ReferenceNorm  = DataType();

      if( this->Is(LocalFlagType::INCREMENTAL) ){
        CalculateIncrementalNorm(rDofSet,rDx, ReferenceNorm, CorrectionNorm);
      }
      else{
        CalculateReferenceNorm(rDofSet,rDx, ReferenceNorm, CorrectionNorm);
      }

      if(CorrectionNorm != 0)
      {
        ratio = CorrectionNorm/ReferenceNorm;
      }
      
      if( ratio == 0 && int(CorrectionNorm-int(ReferenceNorm))==0 )
      {
        ratio = 1.0;
      }
      
      const double Dx_size = TSparseSpace::Size(rDx);

      const double absolute_norm = (CorrectionNorm/std::sqrt(Dx_size));

      if (rModelPart.GetCommunicator().MyPID() == 0)
      {
        if(this->GetEchoLevel() >= 1)
        {
          std::cout << "VARIABLE CRITERION :: Ratio = " << ratio << "; Norm = " << absolute_norm << std::endl;
        }
      }

      rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
      rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

      if ( ratio <= mRatioTolerance  ||  absolute_norm < mAlwaysConvergedNorm )
      {
        if (rModelPart.GetCommunicator().MyPID() == 0 )
        {
          if( this->GetEchoLevel() >= 1 )
          {
            std::cout << "Convergence is achieved" << std::endl;
          }
        }
        return true;
      }
      else
      {
        if( this->Is(LocalFlagType::INCREMENTAL) && ratio == 1.0  && ReferenceNorm <= mRatioTolerance * 1e-2)
        {
          if (rModelPart.GetCommunicator().MyPID() == 0)
          {
            if (this->GetEchoLevel() >= 1)
            {
              std::cout << "Convergence is achieved : - no movement - " << std::endl;
            }
          }		  
          return true;
        }
        else
        {
          return false;
        }       
      }
    }
    else //in this case all the dofs are imposed!
    {
      return true;
    }
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

  DataType mAlwaysConvergedNorm;

  ///@}
  ///@name Private Operators
  ///@{

  void CalculateReferenceNorm(DofsArrayType& rDofSet, const SystemVectorType& rDx, DataType& rReferenceNorm, DataType& rCorrectionNorm)
  {

    rCorrectionNorm = TSparseSpace::TwoNorm(rDx);
    
    DataType temp;

    for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
    {
      if(i_dof->IsFree())
      {
        temp = i_dof->GetSolutionStepValue();
        rReferenceNorm += temp*temp;
      }
    }

    rReferenceNorm = sqrt(rReferenceNorm+1e-20); //to avoid 0
  }


  void CalculateIncrementalNorm(DofsArrayType& rDofSet, const SystemVectorType& rDx, DataType& rReferenceNorm, DataType& rCorrectionNorm)
  {
    DataType temp;
    
    for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
    {
      if(i_dof->IsFree())
      {
        //here we do: d_n+1^it-d_n
        temp = rDx[i_dof->EquationId()];
        rCorrectionNorm +=  temp*temp;
        temp = (i_dof->GetSolutionStepValue()-i_dof->GetSolutionStepValue(1));
        rReferenceNorm += temp*temp;
      }
    }

    rCorrectionNorm = sqrt(rCorrectionNorm);
    rReferenceNorm = sqrt(rReferenceNorm+1e-20); //to avoid 0
  }

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

}; // Class VariableCriterion

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_VARIABLE_CRITERION_H_INCLUDED  defined 

