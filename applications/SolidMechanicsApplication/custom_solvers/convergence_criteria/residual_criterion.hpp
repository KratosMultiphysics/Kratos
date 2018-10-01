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

  typedef std::vector<VariableComponent<VectorComponentType> > ComponentVariableVector;

  /// Pointer definition of ResidualCriterion
  KRATOS_CLASS_POINTER_DEFINITION(ResidualCriterion);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Constructor.
  ResidualCriterion(DataType RatioTolerance,
                    DataType AbsoluteTolerance)
      : BaseType(), mRatioTolerance(RatioTolerance), mAbsoluteTolerance(AbsoluteTolerance)
  {
    mpScalarVariable = nullptr;
    mpVectorVariable = nullptr;
    this->Set(LocalFlagType::SUPPLIED_DOF, false);
  }

  /// Constructor.
  ResidualCriterion(const VariableScalarType& rScalarVariable,
                    DataType RatioTolerance,
                    DataType AbsoluteTolerance)
      : BaseType(), mRatioTolerance(RatioTolerance), mAbsoluteTolerance(AbsoluteTolerance)
  {
    mpScalarVariable = &rScalarVariable;
    mpVectorVariable = nullptr;
    this->Set(LocalFlagType::SUPPLIED_DOF, true);
  }

  /// Constructor.
  ResidualCriterion(const VariableVectorType& rVectorVariable,
                    DataType RatioTolerance,
                    DataType AbsoluteTolerance)
      :BaseType(), mRatioTolerance(RatioTolerance), mAbsoluteTolerance(AbsoluteTolerance)
  {
    mpVectorVariable = &rVectorVariable;
    mpScalarVariable = nullptr;
    this->Set(LocalFlagType::SUPPLIED_DOF, true);
  }

  /// Copy constructor.
  ResidualCriterion( ResidualCriterion const& rOther )
      :BaseType(rOther)
      ,mRatioTolerance(rOther.mRatioTolerance)
      ,mAbsoluteTolerance(rOther.mAbsoluteTolerance)
      ,mInitialResidualNorm(rOther.mInitialResidualNorm)
      ,mpScalarVariable(rOther.mpScalarVariable)
      ,mpVectorVariable(rOther.mpVectorVariable)
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

      std::size_t size = 0;

      if( this->IsNot(LocalFlagType::INITIALIZED) )
      {
        size = CalculateResidualNorm(rDofSet,rb,mInitialResidualNorm);
        this->Set(LocalFlagType::INITIALIZED,true);
      }

      DataType ratio;
      DataType CurrentResidualNorm;
      size = CalculateResidualNorm(rDofSet,rb,CurrentResidualNorm);

      if( size == 0 )
        KRATOS_ERROR << "Dofs vector has size: " << size << std::endl;

      if(mInitialResidualNorm == 0.00)
      {
        ratio = 0.00;
      }
      else
      {
        ratio = CurrentResidualNorm/mInitialResidualNorm;
      }

      const DataType absolute_norm = (CurrentResidualNorm/static_cast<DataType>(size));

      if (rModelPart.GetCommunicator().MyPID() == 0)
      {
        if (this->GetEchoLevel() >= 1)
        {
          std::cout << "RESIDUAL (" << GetDofName() << ") ["<<rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER]<<"] :: Ratio = "<< ratio  << "; Norm = " << absolute_norm <<std::endl;
        }
      }

      rModelPart.GetProcessInfo()[CONVERGENCE_RATIO] = ratio;
      rModelPart.GetProcessInfo()[RESIDUAL_NORM] = absolute_norm;

      if (ratio <= mRatioTolerance || absolute_norm < mAbsoluteTolerance)
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

  DataType mAbsoluteTolerance;

  DataType mInitialResidualNorm;

  VariableScalarPointer  mpScalarVariable;

  VariableVectorPointer  mpVectorVariable;

  ///@}
  ///@name Private Operators
  ///@{

  ///@}
  ///@name Private Operations
  ///@{

  std::size_t CalculateResidualNorm(DofsArrayType& rDofSet, const SystemVectorType& rb, DataType& rResidualNorm)
  {
    KRATOS_TRY

    std::size_t size = 0;
    rResidualNorm = 0.0;

    if( this->Is(CriterionLocalFlags::SUPPLIED_DOF) ){

      if( mpVectorVariable != nullptr ){

        DataType temp;
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
          if(i_dof->IsFree())
          {
            if(CheckVectorDof(i_dof))
            {
              temp = rb[i_dof->EquationId()];
              rResidualNorm +=  temp*temp;
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
              temp = rb[i_dof->EquationId()];
              rResidualNorm += temp*temp;
              ++size;
            }
          }
        }

      }
      else{
        KRATOS_ERROR << " No variable dof supplied to the convergence criterion " << std::endl;
      }

      rResidualNorm = std::sqrt(rResidualNorm);
    }
    else{
      rResidualNorm = TSparseSpace::TwoNorm(rb);
      size = TSparseSpace::Size(rb);
    }

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
