//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_NEWMARK_STEP_METHOD_H_INCLUDED)
#define  KRATOS_NEWMARK_STEP_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/time_integration_methods/newmark_method.hpp"

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


  /// Short class definition.
  /** Detail class definition.
   * This class performs predict and update of dofs variables, their time derivatives and time integrals
   */
  template<class TVariableType, class TValueType>
  class NewmarkStepMethod : public NewmarkMethod<TVariableType,TValueType>
  {
  public:

    ///@name Type Definitions
    ///@{

    /// BaseType
    typedef TimeIntegrationMethod<TVariableType,TValueType>  BaseType;

    /// BasePointerType
    typedef typename BaseType::Pointer                BasePointerType;

    /// NodeType
    typedef typename BaseType::NodeType                      NodeType;

    /// KratosVariable or KratosVariableComponent
    typedef typename BaseType::VariablePointer        VariablePointer;

    /// DerivedType
    typedef NewmarkMethod<TVariableType,TValueType>       DerivedType;


    KRATOS_CLASS_POINTER_DEFINITION( NewmarkStepMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    NewmarkStepMethod() : DerivedType()
    {
      mpStepVariable = nullptr;
    }

    /// Constructor.
    NewmarkStepMethod(const TVariableType& rVariable) : DerivedType(rVariable)
    {
      mpStepVariable = nullptr;
    }

    /// Constructor.
    NewmarkStepMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative)
    {
      mpStepVariable = nullptr;
    }

    /// Constructor.
    NewmarkStepMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable)
    {
      mpStepVariable = nullptr;
    }

    /// Copy Constructor.
    NewmarkStepMethod(NewmarkStepMethod& rOther)
      :DerivedType(rOther)
      ,mpStepVariable(rOther.mpStepVariable)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new NewmarkStepMethod(*this) );
    }

    /// Destructor.
    ~NewmarkStepMethod() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // has step variable
    bool HasStepVariable() override
    {
      return true;
    }

    // set step variable (step variable)
    void SetStepVariable(const TVariableType& rStepVariable) override
    {
      mpStepVariable = &rStepVariable;
    }


    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * @return 0 all ok
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override
    {
      KRATOS_TRY

      // Perform base integration method checks
      int ErrorCode = 0;
      ErrorCode = BaseType::Check(rCurrentProcessInfo);

      if( this->mpStepVariable == nullptr ){
        KRATOS_ERROR << " time integration method Step Variable not set " <<std::endl;
      }
      else{
        KRATOS_CHECK_VARIABLE_KEY((*this->mpStepVariable));
      }

      return ErrorCode;

      KRATOS_CATCH("")
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NewmarkStepMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NewmarkStepMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NewmarkStepMethod Data";
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

    // method variables
    VariablePointer mpStepVariable;

    ///@}
    ///@name Protected Operators
    ///@{

    void AssignFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from first derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      // newmark consistent
      CurrentVariable = PreviousVariable + this->mNewmark.delta_time * ( ( 1.0 - (this->mNewmark.beta/this->mNewmark.gamma) ) * PreviousFirstDerivative + (this->mNewmark.beta/this->mNewmark.gamma) * CurrentFirstDerivative +  this->mNewmark.delta_time *  0.5 * ( 1.0 - 2.0 * (this->mNewmark.beta/this->mNewmark.gamma) ) * PreviousSecondDerivative );
      // uniform accelerated movement
      //CurrentVariable = PreviousVariable + 0.5 * mNewmark.delta_time * (PreviousFirstDerivative + CurrentFirstDerivative) + 0.5 * mNewmark.delta_time * mNewmark.delta_time * PreviousSecondDerivative;

      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      // variable prediction :: uniform accelerated movement **
      // CurrentVariable -= this->mNewmark.delta_time * PreviousFirstDerivative;

      CurrentFirstDerivative   = PreviousFirstDerivative;
      CurrentSecondDerivative -= CurrentSecondDerivative;

      KRATOS_CATCH( "" )
    }

    void AssignFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from second derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      CurrentVariable = PreviousVariable + this->mNewmark.delta_time * PreviousFirstDerivative + this->mNewmark.delta_time * this->mNewmark.delta_time * ( 0.5 * ( 1.0 - 2.0 * this->mNewmark.beta ) * PreviousSecondDerivative + this->mNewmark.beta * CurrentSecondDerivative );

      // variable prediction :: uniform accelerated movement **
      // CurrentVariable -= this->mNewmark.delta_time * PreviousFirstDerivative + 0.5 * this->mNewmark.delta_time * this->mNewmark.delta_time * PreviousSecondDerivative;

      CurrentSecondDerivative = PreviousSecondDerivative;

      KRATOS_CATCH( "" )
    }

    ///@}
    ///@name Protected Operations
    ///@{


    void PredictFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      this->PredictStepVariable(rNode);
      this->PredictFirstDerivative(rNode);
      this->PredictSecondDerivative(rNode);
      this->PredictVariable(rNode);

      // const TValueType& CurrentVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,     0);
      // const TValueType& CurrentStepVariable       = rNode.FastGetSolutionStepValue(*this->mpStepVariable, 0);
      // const TValueType& CurrentFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
      // const TValueType& CurrentSecondDerivative   = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      // std::cout<<*this->mpVariable<<" Predict Node["<<rNode.Id()<<"]"<<CurrentVariable<<" "<<CurrentStepVariable<<" "<<CurrentFirstDerivative<<" "<<CurrentSecondDerivative<<std::endl;

      KRATOS_CATCH( "" )
    }

    virtual void PredictStepVariable(NodeType& rNode)
    {
      KRATOS_TRY

      // predict step variable from previous and current values
      // TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
      // const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      // const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // CurrentStepVariable  = CurrentVariable-PreviousVariable;

      KRATOS_CATCH( "" )
    }

    void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable            = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& PreviousVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // update variable previous iteration instead of previous step
      PreviousVariable = CurrentVariable;

      KRATOS_CATCH( "" )
    }

    void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentFirstDerivative = this->mNewmark.c1 * CurrentStepVariable - this->mNewmark.c4 * CurrentFirstDerivative - this->mNewmark.c5 * CurrentSecondDerivative;

      KRATOS_CATCH( "" )
    }

    void PredictSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      CurrentSecondDerivative  = ( 1.0 / (this->mNewmark.gamma * this->mNewmark.delta_time) ) * ( CurrentFirstDerivative - PreviousFirstDerivative - ( 1.0 - this->mNewmark.gamma ) * this->mNewmark.delta_time * CurrentSecondDerivative );

      KRATOS_CATCH( "" )
    }


    void UpdateFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      this->UpdateStepVariable(rNode);
      this->UpdateFirstDerivative(rNode);
      this->UpdateSecondDerivative(rNode);
      this->UpdateVariable(rNode);

      // const TValueType& CurrentVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,     0);
      // const TValueType& CurrentStepVariable       = rNode.FastGetSolutionStepValue(*this->mpStepVariable, 0);
      // const TValueType& CurrentFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
      // const TValueType& CurrentSecondDerivative   = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      // std::cout<<*this->mpVariable<<" Update Node["<<rNode.Id()<<"]"<<CurrentVariable<<" "<<CurrentStepVariable<<" "<<CurrentFirstDerivative<<" "<<CurrentSecondDerivative<<std::endl;

      KRATOS_CATCH( "" )
    }


    virtual void UpdateStepVariable(NodeType& rNode)
    {
      KRATOS_TRY

      // predict step variable from previous and current values
      TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      CurrentStepVariable += CurrentVariable-PreviousVariable;

      KRATOS_CATCH( "" )
    }

    void UpdateVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& PreviousVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // update variable previous iteration instead of previous step
      PreviousVariable = CurrentVariable;

      KRATOS_CATCH( "" )
    }

    void UpdateFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      CurrentFirstDerivative = (this->mNewmark.c1 * (CurrentStepVariable) - this->mNewmark.c4 * PreviousFirstDerivative - this->mNewmark.c5 * PreviousSecondDerivative);

      KRATOS_CATCH( "" )
    }

    void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentStepVariable      = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      CurrentSecondDerivative = (this->mNewmark.c0 * (CurrentStepVariable) - this->mNewmark.c2 * PreviousFirstDerivative - this->mNewmark.c3 * PreviousSecondDerivative);

      KRATOS_CATCH( "" )
    }

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
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
      // rSerializer.save("StepVariable", mpStepVariable);
    };

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      // rSerializer.load("StepVariable", mpStepVariable);
    };

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NewmarkStepMethod

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, NewmarkStepMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const NewmarkStepMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEWMARK_STEP_METHOD_H_INCLUDED defined
