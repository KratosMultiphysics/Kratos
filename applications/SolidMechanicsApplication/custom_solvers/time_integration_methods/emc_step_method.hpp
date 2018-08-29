//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EMC_STEP_METHOD_H_INCLUDED)
#define  KRATOS_EMC_STEP_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/time_integration_methods/time_integration_method.hpp"

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
  class EmcStepMethod : public TimeIntegrationMethod<TVariableType,TValueType>
  {
  protected:

    struct EmcParameters
    {
      double alpha;
      double delta_time;

      //system constants
      double c0;
      double c1;

      void SetParameters(const double& ralpha,
			 const double& rdelta_time)
      {
	alpha = ralpha;

	delta_time = rdelta_time;

	c0 = ( 2.0 / delta_time );
        c1 = ( 1.0 / delta_time );
      }


    private:

      friend class Serializer;

      void save(Serializer& rSerializer) const
      {
	rSerializer.save("alpha", alpha);
	rSerializer.save("delta_time", delta_time);
	rSerializer.save("c0", c0);
	rSerializer.save("c1", c1);
      };

      void load(Serializer& rSerializer)
      {
	rSerializer.load("alpha", alpha);
	rSerializer.load("delta_time", delta_time);
	rSerializer.load("c0", c0);
	rSerializer.load("c1", c1);
      };

    };

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

    KRATOS_CLASS_POINTER_DEFINITION( EmcStepMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    EmcStepMethod() : BaseType()
    {
      mpStepVariable = nullptr;
    }

    /// Constructor.
    EmcStepMethod(const TVariableType& rVariable) : BaseType(rVariable)
    {
      mpStepVariable = nullptr;
    }

    /// Constructor.
    EmcStepMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : BaseType(rVariable,rFirstDerivative,rSecondDerivative)
    {
      mpStepVariable = nullptr;
    }

    /// Constructor.
    EmcStepMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : BaseType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable)
    {
      mpStepVariable = nullptr;
    }

    /// Copy Constructor.
    EmcStepMethod(EmcStepMethod& rOther)
      :BaseType(rOther)
      ,mEmc(rOther.mEmc)
      ,mpStepVariable(rOther.mpStepVariable)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new EmcStepMethod(*this) );
    }

    /// Destructor.
    ~EmcStepMethod() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    //calculate parameters (to call it once with the original input parameters)
    void CalculateParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY


     double alpha = 0.5;
     if (rCurrentProcessInfo.Has(EQUILIBRIUM_POINT))
       {
	 alpha = rCurrentProcessInfo[EQUILIBRIUM_POINT];
       }

     rCurrentProcessInfo[EQUILIBRIUM_POINT] = alpha;

     this->SetParameters(rCurrentProcessInfo);

     KRATOS_CATCH( "" )
    }


    // set parameters (do not calculate parameters here, only read them)
    void SetParameters(const ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     double delta_time = rCurrentProcessInfo[DELTA_TIME];

     if (delta_time < 1.0e-24)
        {
	  KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Method DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part " << std::endl;
        }

     double alpha = 0.5;
     if (rCurrentProcessInfo.Has(EQUILIBRIUM_POINT))
       {
	 alpha = rCurrentProcessInfo[EQUILIBRIUM_POINT];
       }

     mEmc.SetParameters(alpha,delta_time);

     KRATOS_CATCH( "" )
    }

    // set parameters to process info
    void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     rCurrentProcessInfo[EQUILIBRIUM_POINT]  = this->mEmc.alpha;

     KRATOS_CATCH( "" )
    }

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
        buffer << "EmcStepMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmcStepMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "EmcStepMethod Data";
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

    // method parameters
    EmcParameters   mEmc;

    // method variables
    VariablePointer mpStepVariable;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    void AssignFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from variable


      KRATOS_CATCH( "" )
    }


    void AssignFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from first derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      TValueType& PreviousFirstDerivative        = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      PreviousFirstDerivative = CurrentFirstDerivative;

      CurrentVariable = PreviousVariable + (CurrentFirstDerivative+PreviousFirstDerivative) * (1.0/this->mEmc.c0);

      KRATOS_CATCH( "" )
    }

    void AssignFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from second derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentVariable = PreviousVariable + CurrentSecondDerivative * (1.0/(this->mEmc.c0*this->mEmc.c1));

      KRATOS_CATCH( "" )
    }

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
      TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      CurrentStepVariable  = CurrentVariable-PreviousVariable;

      KRATOS_CATCH( "" )
    }

    void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& PreviousVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // update variable previous iteration instead of previous step
      PreviousVariable = CurrentVariable;

      KRATOS_CATCH( "" )
    }

    void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      CurrentFirstDerivative = this->mEmc.c0 * (CurrentVariable-PreviousVariable) - PreviousFirstDerivative;

      KRATOS_CATCH( "" )
    }

    void PredictSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY


      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      CurrentSecondDerivative = this->mEmc.c1 * (CurrentFirstDerivative-PreviousFirstDerivative);


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

      TValueType& CurrentFirstDerivative        = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      const TValueType& CurrentStepVariable     = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);

      CurrentFirstDerivative = this->mEmc.c0 * (CurrentStepVariable) - PreviousFirstDerivative;


      KRATOS_CATCH( "" )
    }

    void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      CurrentSecondDerivative = this->mEmc.c1 * (CurrentFirstDerivative - PreviousFirstDerivative);

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
      rSerializer.save("EmcParameters", mEmc);
      // rSerializer.save("StepVariable", mpStepVariable);
    };

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      rSerializer.load("EmcParameters", mEmc);
      // rSerializer.load("StepVariable", mpStepVariable);
    };

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class EmcStepMethod

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, EmcStepMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const EmcStepMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_EMC_STEP_METHOD_H_INCLUDED defined
