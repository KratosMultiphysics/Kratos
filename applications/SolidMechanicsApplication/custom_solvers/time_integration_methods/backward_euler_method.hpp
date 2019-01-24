//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               April 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BACKWARD_EULER_METHOD_H_INCLUDED)
#define  KRATOS_BACKWARD_EULER_METHOD_H_INCLUDED

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
  class BackwardEulerMethod : public TimeIntegrationMethod<TVariableType,TValueType>
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

    KRATOS_CLASS_POINTER_DEFINITION( BackwardEulerMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    BackwardEulerMethod() : BaseType() {}

    /// Constructor.
    BackwardEulerMethod(const TVariableType& rVariable) : BaseType(rVariable) {}

    /// Constructor.
    BackwardEulerMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : BaseType(rVariable,rFirstDerivative,rSecondDerivative) {}

    /// Constructor.
    BackwardEulerMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : BaseType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable) {}

    /// Copy Constructor.
    BackwardEulerMethod(BackwardEulerMethod& rOther)
      :BaseType(rOther)
      ,mDeltaTime(rOther.mDeltaTime)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new BackwardEulerMethod(*this) );
    }

    /// Destructor.
    ~BackwardEulerMethod() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // set parameters (do not calculate parameters here, only read them)
    void SetParameters(const ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

     if (delta_time < 1.0e-24)
        {
	  KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Method DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part " << std::endl;
        }

     mDeltaTime = delta_time;

     KRATOS_CATCH( "" )
    }

    // set parameters to process info
    void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY


     KRATOS_CATCH( "" )
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

      // Check that all required variables have been registered
      if( this->mpFirstDerivative == nullptr ){
        KRATOS_ERROR << " time integration method FirstDerivative not set " <<std::endl;
      }
      else{
        KRATOS_CHECK_VARIABLE_KEY((*this->mpFirstDerivative));
      }

      if( this->mpSecondDerivative == nullptr ){
        KRATOS_ERROR << " time integration method SecondDerivative not set " <<std::endl;
      }
      else{
        KRATOS_CHECK_VARIABLE_KEY((*this->mpSecondDerivative));
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
        buffer << "BackwardEulerMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BackwardEulerMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "BackwardEulerMethod Data";
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

    double  mDeltaTime;

    ///@}
    ///@name Protected Operators
    ///@{

    void AssignFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from variable
      TValueType& CurrentFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentFirstDerivative  -= CurrentFirstDerivative;
      CurrentSecondDerivative -= CurrentSecondDerivative;


      KRATOS_CATCH( "" )
    }

    void AssignFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from first derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& OldPreviousVariable      = rNode.FastGetSolutionStepValue(*this->mpVariable,         2);

      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      // backward euler consistent
      CurrentVariable = (1.0/3.0) * ( 4.0 * PreviousVariable - OldPreviousVariable + 2.0 * mDeltaTime * CurrentFirstDerivative);

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

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

      const TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      // backward euler consistent
      CurrentVariable = PreviousVariable + mDeltaTime * ( PreviousFirstDerivative + 0.5 * mDeltaTime * CurrentSecondDerivative );

      KRATOS_CATCH( "" )
    }

    void PredictFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      this->PredictVariable(rNode);
      this->PredictFirstDerivative(rNode);
      this->PredictSecondDerivative(rNode);

      KRATOS_CATCH( "" )
    }

    void PredictFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      this->PredictVariable(rNode);
      this->PredictFirstDerivative(rNode);
      this->PredictSecondDerivative(rNode);

      KRATOS_CATCH( "" )
    }

    void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& PreviousFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      // variable prediction
      CurrentVariable = PreviousVariable + mDeltaTime * PreviousFirstDerivative;


      KRATOS_CATCH( "" )
    }

    void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // TValueType& CurrentFirstDerivative           = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      // const TValueType& PreviousFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      // CurrentFirstDerivative = PreviousFirstDerivative;

      const TValueType& CurrentVariable            = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& OldPreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         2);

      TValueType& CurrentFirstDerivative           = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      CurrentFirstDerivative = (1.0/mDeltaTime) * ( 3.0 * 0.5 * CurrentVariable - 2.0 * PreviousVariable + 0.5 * OldPreviousVariable);

      KRATOS_CATCH( "" )
    }

    void PredictSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // TValueType& CurrentSecondDerivative          = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
      // const TValueType& PreviousSecondDerivative   = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      // CurrentSecondDerivative = PreviousSecondDerivative;

      const TValueType& CurrentVariable            = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& PreviousFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      TValueType& CurrentSecondDerivative          = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentSecondDerivative = (2.0/ (mDeltaTime*mDeltaTime)) * (CurrentVariable - PreviousVariable - mDeltaTime * PreviousFirstDerivative);

      KRATOS_CATCH( "" )
    }

    void UpdateFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable            = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& OldPreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         2);

      TValueType& CurrentFirstDerivative           = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      TValueType& CurrentSecondDerivative          = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentSecondDerivative = (2.0/ (mDeltaTime*mDeltaTime)) * (CurrentVariable - PreviousVariable - mDeltaTime * PreviousFirstDerivative);

      CurrentFirstDerivative = (1.0/mDeltaTime) * ( (3.0/2.0) * CurrentVariable - 2.0 * PreviousVariable + 0.5 * OldPreviousVariable);

      KRATOS_CATCH( "" )
    }


    void UpdateFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY


      TValueType& CurrentVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& OldPreviousVariable     = rNode.FastGetSolutionStepValue(*this->mpVariable,         2);

      const TValueType& CurrentFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      TValueType& CurrentSecondDerivative       = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentVariable = (2.0 / 3.0) * ( mDeltaTime * CurrentFirstDerivative + 2.0 * PreviousVariable - 0.5 * OldPreviousVariable);
      CurrentSecondDerivative = (2.0/ (mDeltaTime*mDeltaTime)) * (CurrentVariable - PreviousVariable - mDeltaTime * PreviousFirstDerivative);


      KRATOS_CATCH( "" )
    }

    void UpdateFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable                  = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& OldPreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         2);

      TValueType& CurrentFirstDerivative           = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      const TValueType& CurrentSecondDerivative    = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentVariable = PreviousVariable + mDeltaTime * PreviousFirstDerivative + 0.5 * mDeltaTime * mDeltaTime * CurrentSecondDerivative;
      CurrentFirstDerivative = (1.0/mDeltaTime) * ( 3.0 * 0.5 * CurrentVariable - 2.0 * PreviousVariable + 0.5 * OldPreviousVariable);

      KRATOS_CATCH( "" )
    }

    void UpdateVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      KRATOS_CATCH( "" )
    }


    void UpdateFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable            = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& OldPreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         2);

      TValueType& CurrentFirstDerivative           = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      CurrentFirstDerivative = (1.0/mDeltaTime) * ( 3.0 * 0.5 * CurrentVariable - 2.0 * PreviousVariable + 0.5 * OldPreviousVariable);

      KRATOS_CATCH( "" )
    }

    void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable            = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& PreviousFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      TValueType& CurrentSecondDerivative          = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentSecondDerivative = (2.0/ (mDeltaTime*mDeltaTime)) * (CurrentVariable - PreviousVariable - mDeltaTime * PreviousFirstDerivative);

      KRATOS_CATCH( "" )
    }

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    // get parameters
    double& GetFirstDerivativeInertialParameter(double& rParameter) override
    {
      rParameter = 3.0 / (2.0 * mDeltaTime);
      return rParameter;
    }

    double& GetSecondDerivativeInertialParameter(double& rParameter) override
    {
      rParameter = 2.0 / (mDeltaTime*mDeltaTime);
      return rParameter;
    }

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
      rSerializer.save("DeltaTime", mDeltaTime);
    };

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      rSerializer.load("DeltaTime", mDeltaTime);
    };

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class BackwardEulerMethod

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, BackwardEulerMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const BackwardEulerMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BACKWARD_EULER_METHOD_H_INCLUDED defined
