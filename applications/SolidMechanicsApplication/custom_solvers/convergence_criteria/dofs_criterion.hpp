//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DOFS_CRITERION_H_INCLUDED )
#define  KRATOS_DOFS_CRITERION_H_INCLUDED

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


/** @class DofsCriterion
 * @brief This convergence criteria checks the variable dofs
 */

template<class TSparseSpace, class TDenseSpace>
class DofsCriterion : public ConvergenceCriterion<TSparseSpace,TDenseSpace>
{
 public:

  ///@name Type Definitions
  ///@{

  typedef ConvergenceCriterion< TSparseSpace, TDenseSpace >     BaseType;
  typedef typename BaseType::LocalFlagType                 LocalFlagType;
  typedef typename BaseType::DataType                           DataType;
  typedef typename BaseType::DofsArrayType                 DofsArrayType;
  typedef typename BaseType::SystemMatrixType           SystemMatrixType;
  typedef typename BaseType::SystemVectorType           SystemVectorType;

  typedef array_1d<double,3>                                  VectorType;
  typedef VectorComponentAdaptor<VectorType>         VectorComponentType;
  typedef Variable<VectorType>                        VariableVectorType;
  typedef Variable<double>                            VariableScalarType;
  typedef const VariableVectorType*                VariableVectorPointer;
  typedef const VariableScalarType*                VariableScalarPointer;

  /// Pointer definition of DofsCriterion
  KRATOS_CLASS_POINTER_DEFINITION( DofsCriterion );

  ///@}
  ///@name Life Cycle
  ///@{

  /// Constructor.
  DofsCriterion(DataType RatioTolerance,
                DataType AbsoluteTolerance)
      : BaseType(), mRatioTolerance(RatioTolerance), mAbsoluteTolerance(AbsoluteTolerance)
  {
    mpScalarVariable = nullptr;
    mpVectorVariable = nullptr;
    this->Set(LocalFlagType::SUPPLIED_DOF, false);
  }

  /// Constructor.
  DofsCriterion(const VariableScalarType& rScalarVariable,
                DataType RatioTolerance,
                DataType AbsoluteTolerance)
      : BaseType(), mRatioTolerance(RatioTolerance), mAbsoluteTolerance(AbsoluteTolerance)
  {
    mpScalarVariable = &rScalarVariable;
    mpVectorVariable = nullptr;
    this->Set(LocalFlagType::SUPPLIED_DOF, true);
  }

  /// Constructor.
  DofsCriterion(const VariableVectorType& rVectorVariable,
                DataType RatioTolerance,
                DataType AbsoluteTolerance)
      : BaseType(), mRatioTolerance(RatioTolerance), mAbsoluteTolerance(AbsoluteTolerance)
  {
    mpVectorVariable = &rVectorVariable;
    mpScalarVariable = nullptr;
    this->Set(LocalFlagType::SUPPLIED_DOF, true);
  }

  /// Copy constructor.
  DofsCriterion(DofsCriterion const& rOther)
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAbsoluteTolerance(rOther.mAbsoluteTolerance)
      ,mpScalarVariable(rOther.mpScalarVariable)
      ,mpVectorVariable(rOther.mpVectorVariable)
  {
  }

  /// Destructor.
  ~DofsCriterion() override {}


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

      std::size_t size = 0;

      if( this->Is(LocalFlagType::INCREMENTAL) ){
        size = CalculateIncrementalNorm(rDofSet,rDx, ReferenceNorm, CorrectionNorm);
      }
      else{
        size = CalculateReferenceNorm(rDofSet,rDx, ReferenceNorm, CorrectionNorm);
      }

      if( size == 0 ){
        KRATOS_WARNING("") << GetDofName() <<" Dofs vector has size: " << size << std::endl;
        return true;
      }

      if(CorrectionNorm != 0)
      {
        ratio = CorrectionNorm/ReferenceNorm;
      }

      if( ratio == 0 && int(CorrectionNorm-int(ReferenceNorm))==0 )
      {
        ratio = 1.0;
      }

      const DataType absolute_norm = (CorrectionNorm/static_cast<DataType>(size));

      if (rModelPart.GetCommunicator().MyPID() == 0)
      {
        if(this->GetEchoLevel() >= 1)
        {
          std::cout << "DOF (" << GetDofName() << ") ["<<rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER]<<"] :: Ratio = " << ratio << "; Norm = " << absolute_norm << std::endl;
          std::cout << " CorrectionNorm = " << CorrectionNorm << "; ReferenceNorm = " << ReferenceNorm << std::endl;
        }
      }

      rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
      rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

      if ( ratio <= mRatioTolerance  ||  absolute_norm < mAbsoluteTolerance )
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

  DataType mAbsoluteTolerance;

  VariableScalarPointer  mpScalarVariable;

  VariableVectorPointer  mpVectorVariable;

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  std::size_t CalculateReferenceNorm(DofsArrayType& rDofSet, const SystemVectorType& rDx, DataType& rReferenceNorm, DataType& rCorrectionNorm)
  {
    KRATOS_TRY

    std::size_t size = 0;
    rCorrectionNorm = 0.0;
    rReferenceNorm  = 0.0;

    if( this->Is(CriterionLocalFlags::SUPPLIED_DOF) ){

      if( mpVectorVariable != nullptr ){

        DataType temp;
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
          if(i_dof->IsFree())
          {
            if(CheckVectorDof(i_dof))
            {
              temp = rDx[i_dof->EquationId()];
              rCorrectionNorm +=  temp*temp;
              temp = i_dof->GetSolutionStepValue();
              rReferenceNorm += temp*temp;
              ++size;
            }
          }
        }
      }
      else if( mpScalarVariable != nullptr ){

        DataType temp;

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
          if(i_dof->IsFree())
          {
            if(i_dof->GetVariable() == *mpScalarVariable)
            {
              temp = rDx[i_dof->EquationId()];
              rCorrectionNorm +=  temp*temp;
              temp = i_dof->GetSolutionStepValue();
              rReferenceNorm += temp*temp;
              ++size;
            }
          }
        }

      }
      else{
        KRATOS_ERROR << " No variable dof supplied to the convergence criterion " << std::endl;
      }

      rCorrectionNorm = std::sqrt(rCorrectionNorm);
    }
    else{

      rCorrectionNorm = TSparseSpace::TwoNorm(rDx);
      size = TSparseSpace::Size(rDx);

      DataType temp;

      for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
      {
        if(i_dof->IsFree())
        {
          temp = i_dof->GetSolutionStepValue();
          rReferenceNorm += temp*temp;
        }
      }

    }

    rReferenceNorm = std::sqrt(rReferenceNorm+1e-20); //to avoid 0

    return size;

    KRATOS_CATCH("")
  }


  std::size_t CalculateIncrementalNorm(DofsArrayType& rDofSet, const SystemVectorType& rDx, DataType& rReferenceNorm, DataType& rCorrectionNorm)
  {
    KRATOS_TRY

    std::size_t size = 0;
    rCorrectionNorm = 0.0;
    rReferenceNorm  = 0.0;

    if( this->Is(CriterionLocalFlags::SUPPLIED_DOF) ){

      if( mpVectorVariable != nullptr ){

        DataType temp;

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
          if(i_dof->IsFree())
          {
            if(CheckVectorDof(i_dof))
            {
              temp = rDx[i_dof->EquationId()];
              rCorrectionNorm +=  temp*temp;
              temp = (i_dof->GetSolutionStepValue()-i_dof->GetSolutionStepValue(1));
              rReferenceNorm += temp*temp;
              ++size;
            }
          }
        }
      }
      else if( mpScalarVariable != nullptr ){

        DataType temp;

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
          if(i_dof->IsFree())
          {
            if(i_dof->GetVariable() == *mpScalarVariable)
            {
              temp = rDx[i_dof->EquationId()];
              rCorrectionNorm +=  temp*temp;
              temp = (i_dof->GetSolutionStepValue()-i_dof->GetSolutionStepValue(1));
              rReferenceNorm += temp*temp;
              ++size;
            }
          }
        }

      }
      else{
        KRATOS_ERROR << " No variable dof supplied to the convergence criterion " << std::endl;
      }
    }
    else{

      DataType temp;

      for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
      {
        if(i_dof->IsFree())
        {
          temp = rDx[i_dof->EquationId()];
          rCorrectionNorm +=  temp*temp;
          temp = (i_dof->GetSolutionStepValue()-i_dof->GetSolutionStepValue(1));
          rReferenceNorm += temp*temp;
          ++size;
        }
      }

    }

    rCorrectionNorm = std::sqrt(rCorrectionNorm);
    rReferenceNorm  = std::sqrt(rReferenceNorm+1e-20); //to avoid 0

    return size;

    KRATOS_CATCH("")
  }

  bool CheckVectorDof(typename DofsArrayType::iterator& rDofIter)
  {
    KRATOS_TRY

    const std::string& variable_name = mpVectorVariable->Name();
    const VariableComponent<VectorComponentType>& var_x = KratosComponents<VariableComponent<VectorComponentType> >::Get(variable_name+"_X");
    const VariableComponent<VectorComponentType>& var_y = KratosComponents<VariableComponent<VectorComponentType> >::Get(variable_name+"_Y");
    const VariableComponent<VectorComponentType>& var_z = KratosComponents<VariableComponent<VectorComponentType> >::Get(variable_name+"_Z");

    if( rDofIter->GetVariable() == var_x || rDofIter->GetVariable() == var_y || rDofIter->GetVariable() == var_z )
      return true;
    else
      return false;

    KRATOS_CATCH("")
  }

  std::string GetDofName()
  {
    std::string name = "DOFS";

    if( this->Is(CriterionLocalFlags::SUPPLIED_DOF) ){

      if( mpVectorVariable != nullptr ){
        name =  mpVectorVariable->Name();
      }
      else if( mpScalarVariable != nullptr ){
        name = mpScalarVariable->Name();
      }
    }

    return name;
  }

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

}; // Class DofsCriterion

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DOFS_CRITERION_H_INCLUDED  defined
