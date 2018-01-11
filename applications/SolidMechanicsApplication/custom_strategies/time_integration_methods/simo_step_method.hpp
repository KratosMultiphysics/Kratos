//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SIMO_STEP_METHOD )
#define  KRATOS_SIMO_STEP_METHOD

// System includes

// External includes

// Project includes
#include "custom_strategies/time_integration_methods/bossak_step_method.hpp"

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) SimoStepMethod : public BossakStepMethod<TVariableType,TValueType>
  {   
  public:
 
    ///@name Type Definitions
    ///@{

    /// BaseType
    typedef TimeIntegrationMethod<TVariableType,TValueType>  BaseType;

    /// BaseTypePointer
    typedef typename BaseType::Pointer                BaseTypePointer;
    
    /// NodeType
    typedef typename BaseType::NodeType                      NodeType;
    
    /// KratosVariable or KratosVariableComponent    
    typedef typename BaseType::VariablePointer        VariablePointer;

    /// DerivedType
    typedef BossakStepMethod<TVariableType,TValueType>    DerivedType;

    
    KRATOS_CLASS_POINTER_DEFINITION( SimoStepMethod );

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    SimoStepMethod() : DerivedType() {}

    /// Copy Constructor.
    SimoStepMethod(SimoStepMethod& rOther) : DerivedType(rOther) {}

    /// Clone.
    BaseTypePointer Clone()
    {
      return BaseTypePointer( new SimoStepMethod(*this) );
    }

    /// Destructor.
    ~SimoStepMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    
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

       this->PredictFirstDerivative(rNode);
       this->PredictSecondDerivative(rNode);
       this->PredictStepVariable(rNode);
       this->PredictVariable(rNode);
       
     }       

     // const TValueType& CurrentVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,     0);
     // const TValueType& CurrentStepVariable       = rNode.FastGetSolutionStepValue(*this->mpStepVariable, 0);
     // const TValueType& CurrentFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
     // const TValueType& CurrentSecondDerivative   = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

     // std::cout<<*this->mpVariable<<" Predict Node["<<rNode.Id()<<"]"<<CurrentVariable<<" "<<CurrentStepVariable<<" "<<CurrentFirstDerivative<<" "<<CurrentSecondDerivative<<std::endl;

     
     KRATOS_CATCH( "" )
    }

    
    // update
    virtual void Update(NodeType& rNode) override
    {
     KRATOS_TRY

     this->UpdateFirstDerivative(rNode);
     this->UpdateSecondDerivative(rNode);

     this->UpdateStepVariable(rNode);
     this->UpdateVariable(rNode);
     
     // const TValueType& CurrentVariable           = rNode.FastGetSolutionStepValue(*this->mpVariable,     0);
     // const TValueType& CurrentStepVariable       = rNode.FastGetSolutionStepValue(*this->mpStepVariable, 0);
     // const TValueType& CurrentFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
     // const TValueType& CurrentSecondDerivative   = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

     // std::cout<<*this->mpVariable<<" Update Node["<<rNode.Id()<<"]"<<CurrentVariable<<" "<<CurrentStepVariable<<" "<<CurrentFirstDerivative<<" "<<CurrentSecondDerivative<<std::endl;
     
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
        buffer << "SimoStepMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimoStepMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SimoStepMethod Data";     
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

    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void PredictFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY
	
      // predict variable from variable

      
      KRATOS_CATCH( "" )      
    } 
    
    virtual void PredictFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
	
      // predict variable from first derivative	
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);      

      CurrentVariable  = PreviousVariable + CurrentFirstDerivative  * (1.0/this->mNewmark.c1);
            
      KRATOS_CATCH( "" )      
    }

    virtual void PredictFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
	
      // predict variable from second derivative
      TValueType& CurrentVariable                 = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& CurrentFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);      

      const TValueType& CurrentSecondDerivative   = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      
      CurrentVariable  = PreviousVariable + CurrentFirstDerivative * (1.0/this->mNewmark.c1) + CurrentSecondDerivative * (1.0/this->mNewmark.c0);
 
      
      KRATOS_CATCH( "" )      
    }

    virtual void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& PreviousVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
	
      // update variable previous iteration instead of previous step
      PreviousVariable = CurrentVariable;
      
      KRATOS_CATCH( "" )
    }
    
    virtual void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
	
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
 	          
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      CurrentFirstDerivative = this->mNewmark.c1 * (CurrentVariable-PreviousVariable);
      	
      KRATOS_CATCH( "" )      
    }

    virtual void PredictSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
 
      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      
      CurrentSecondDerivative = (this->mNewmark.c0/this->mNewmark.c1) * (CurrentFirstDerivative-PreviousFirstDerivative);

      
      KRATOS_CATCH( "" )              
    }


    virtual void PredictStepVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict step variable from previous and current values
      TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
      TValueType& PreviousStepVariable           = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     1);
      
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // update step variable previous iteration instead of previous step
      PreviousStepVariable = CurrentStepVariable;
      
      CurrentStepVariable  = CurrentVariable-PreviousVariable;
      
      KRATOS_CATCH( "" )
    }

    virtual void UpdateStepVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict step variable from previous and current values
      TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
      TValueType& PreviousStepVariable           = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     1);
      
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // update step variable previous iteration instead of previous step
      PreviousStepVariable = CurrentStepVariable;
      
      CurrentStepVariable += CurrentVariable-PreviousVariable;
      
      KRATOS_CATCH( "" )
    }
    

    virtual void UpdateVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& PreviousVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
	
      // update variable previous iteration instead of previous step
      PreviousVariable = CurrentVariable;     
	
      KRATOS_CATCH( "" )
    }

    virtual void UpdateFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
	
      const TValueType& CurrentVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative        = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
 	          
      const TValueType& PreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      CurrentFirstDerivative += this->mNewmark.c1 * (CurrentVariable-PreviousVariable);
      
      KRATOS_CATCH( "" )      
    }

    virtual void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
 	          
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      CurrentSecondDerivative += this->mNewmark.c0 * (CurrentVariable-PreviousVariable);
      
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

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
    };

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )
    };
    
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@}
  
  }; // Class SimoStepMethod
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{
  
  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, SimoStepMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const SimoStepMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }
  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_SIMO_STEP_METHOD defined
