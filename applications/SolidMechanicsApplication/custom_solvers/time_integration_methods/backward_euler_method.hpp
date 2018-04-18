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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) BackwardEulerMethod : public TimeIntegrationMethod<TVariableType,TValueType>
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
    virtual ~BackwardEulerMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    
    // set parameters (do not calculate parameters here, only read them)
    virtual void SetParameters(const ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY
       
     double delta_time = rCurrentProcessInfo[DELTA_TIME];

     if (delta_time < 1.0e-24)
        {
	  KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Method DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part " << std::endl;
        }
     

     mDeltaTime = delta_time;
     
     KRATOS_CATCH( "" )
    }     

    // set parameters to process info
    virtual void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY
       
           
     KRATOS_CATCH( "" )
    }     

    
    // get parameters
    virtual double& GetFirstDerivativeParameter(double& rParameter) override
    {
      rParameter = mDeltaTime;
      return rParameter;
    }

    
    // predict
    virtual void Predict(NodeType& rNode) override
    {
     KRATOS_TRY
     
     if( this->mpInputVariable != nullptr ){ 

       if( *this->mpInputVariable == *this->mpVariable ){
	 this->PredictFromVariable(rNode);
       }

       if( *this->mpInputVariable == *this->mpFirstDerivative ){
	 this->PredictFromFirstDerivative(rNode);
       }
       
       if( *this->mpInputVariable == *this->mpSecondDerivative ){
	 this->PredictFromSecondDerivative(rNode);
       }

     }
     else{

       this->PredictVariable(rNode);
       this->PredictFirstDerivative(rNode);
       this->PredictSecondDerivative(rNode);
       
     }
	
     KRATOS_CATCH( "" )
    }

    
    // update
    virtual void Update(NodeType& rNode) override
    {
     KRATOS_TRY

     if( this->mpOutputVariable != nullptr ){
       
       if( *this->mpOutputVariable != *this->mpVariable ){
	 this->UpdateVariable(rNode);
       }

       if( *this->mpOutputVariable != *this->mpFirstDerivative ){
	 this->UpdateFirstDerivative(rNode);
       }
       
       if( *this->mpOutputVariable != *this->mpSecondDerivative ){
	 this->UpdateSecondDerivative(rNode);
       }
     }
     else{
       
       this->UpdateVariable(rNode);
       this->UpdateFirstDerivative(rNode);
       this->UpdateSecondDerivative(rNode);
     }
     
     KRATOS_CATCH( "" )
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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "BackwardEulerMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BackwardEulerMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    virtual void PredictFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from variable
      TValueType& CurrentFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentFirstDerivative  -= CurrentFirstDerivative;
      CurrentSecondDerivative -= CurrentSecondDerivative;
      
      
      KRATOS_CATCH( "" )
    }

    virtual void PredictFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from first derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);              
      // backward euler consistent
      CurrentVariable = PreviousVariable + mDeltaTime * CurrentFirstDerivative;
      
      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
     
      CurrentFirstDerivative   = PreviousFirstDerivative;    
      CurrentSecondDerivative -= CurrentSecondDerivative;

      KRATOS_CATCH( "" )      
    }

    virtual void PredictFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from second derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);     
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);       

      const TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
       
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      // backward euler consistent
      CurrentVariable = PreviousVariable + mDeltaTime * ( PreviousFirstDerivative + mDeltaTime * CurrentSecondDerivative );

       
      KRATOS_CATCH( "" )      
    }


    virtual void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& CurrentFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      
      // variable prediction
      CurrentVariable += mDeltaTime * CurrentFirstDerivative;
      
      KRATOS_CATCH( "" )
    }

    virtual void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // variable prediction
      CurrentFirstDerivative +=  mDeltaTime * CurrentSecondDerivative;

      // uniform accelerated movement
      // CurrentFirstDerivative = (2.0/this->mNewmark.delta_time) * (CurrentVariable-PreviousVariable) - this->mNewmark.delta_time * CurrentSecondDerivative - CurrentFirstDerivative;

      
      KRATOS_CATCH( "" )      
    }

    virtual void PredictSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
      
      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      
      CurrentSecondDerivative = ( 1.0 / (this->mNewmark.gamma * this->mNewmark.delta_time) ) * ( CurrentFirstDerivative - PreviousFirstDerivative - ( 1.0 - this->mNewmark.gamma ) * this->mNewmark.delta_time * CurrentSecondDerivative );
      
      KRATOS_CATCH( "" )              
    }
    
    virtual void UpdateVariable(NodeType& rNode) override
    {
      KRATOS_TRY
      
      KRATOS_CATCH( "" )
    }
    

    virtual void UpdateFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
	
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
 	          
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      CurrentFirstDerivative = (this->mNewmark.c1 * (CurrentVariable-PreviousVariable) - this->mNewmark.c4 * PreviousFirstDerivative - this->mNewmark.c5 * PreviousSecondDerivative);
  
      KRATOS_CATCH( "" )      
    }

    virtual void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
 	          
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);
      
      CurrentSecondDerivative = (this->mNewmark.c0 * (CurrentVariable-PreviousVariable) - this->mNewmark.c2 * PreviousFirstDerivative - this->mNewmark.c3 * PreviousSecondDerivative);
    
      KRATOS_CATCH( "" )              
    }

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
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
      rSerializer.save("NewmarkParameters", mNewmark);
    };

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      rSerializer.load("NewmarkParameters", mNewmark);
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
