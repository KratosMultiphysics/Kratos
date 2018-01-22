//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EMC_STEP_METHOD )
#define  KRATOS_EMC_STEP_METHOD

// System includes

// External includes

// Project includes
#include "custom_strategies/time_integration_methods/time_integration_method.hpp"

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) EmcStepMethod : public TimeIntegrationMethod<TVariableType,TValueType>
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

    /// BaseTypePointer
    typedef typename BaseType::Pointer                BaseTypePointer;
    
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

    /// Copy Constructor.
    EmcStepMethod(EmcStepMethod& rOther)
      :BaseType(rOther)
      ,mEmc(rOther.mEmc)
      ,mpStepVariable(rOther.mpStepVariable)
    {
    }

    /// Clone.
    BaseTypePointer Clone()
    {
      return BaseTypePointer( new EmcStepMethod(*this) );
    }

    /// Destructor.
    ~EmcStepMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

   
    // set parameters
    virtual void SetParameters(const ProcessInfo& rCurrentProcessInfo) override
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
    virtual void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY
          
     KRATOS_CATCH( "" )
    } 
    
    // has step variable
    virtual bool HasStepVariable() override
    {
      return true;
    }
    
    // set step variable (step variable)
    virtual void SetStepVariable(const TVariableType& rStepVariable) override
    {
      mpStepVariable = &rStepVariable;
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
       this->UpdateFirstDerivative(rNode);
       this->UpdateSecondDerivative(rNode);
       this->PredictStepVariable(rNode);      
     }
	
     KRATOS_CATCH( "" )
    }

    
    // update
    virtual void Update(NodeType& rNode) override
    {
     KRATOS_TRY
       
     this->UpdateVariable(rNode);
     this->UpdateFirstDerivative(rNode);
     this->UpdateSecondDerivative(rNode);

     this->UpdateStepVariable(rNode);  

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
        buffer << "EmcStepMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmcStepMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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


    virtual void PredictFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict step variable from variable	
      TValueType& PreviousStepVariable           = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     1);
      
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      PreviousStepVariable = CurrentVariable - PreviousVariable;
		
      KRATOS_CATCH( "" )
    }

    
    virtual void PredictFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
	
      // predict variable from first derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
  
      CurrentVariable = PreviousVariable + (CurrentFirstDerivative-PreviousFirstDerivative) * (1.0/this->mEmc.c0);

      TValueType& PreviousStepVariable           = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     1);
      PreviousStepVariable = CurrentVariable - PreviousVariable;

      
      KRATOS_CATCH( "" )      
    }

    virtual void PredictFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from second derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);
      
      CurrentVariable = PreviousVariable + (CurrentSecondDerivative-PreviousSecondDerivative) * (1.0/this->mEmc.c1);

      TValueType& PreviousStepVariable           = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     1);
      PreviousStepVariable = CurrentVariable - PreviousVariable;
      
      KRATOS_CATCH( "" )      
    }


    virtual void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

	
      KRATOS_CATCH( "" )
    }

    
    virtual void PredictStepVariable(NodeType& rNode)
    {
      KRATOS_TRY

      this->UpdateStepVariable(rNode);
	
      KRATOS_CATCH( "" )
    }

   
    virtual void UpdateStepVariable(NodeType& rNode)
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

    virtual void UpdateVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      
      KRATOS_CATCH( "" )
    }

    virtual void UpdateFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY
	
      TValueType& CurrentFirstDerivative        = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      
      
      const TValueType& CurrentVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
 	          
      const TValueType& PreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      CurrentFirstDerivative = this->mEmc.c0 * (CurrentVariable-PreviousVariable) - PreviousFirstDerivative;
      
      KRATOS_CATCH( "" )      
    }

    virtual void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& CurrentFirstDerivative        = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
                  
      CurrentSecondDerivative = this->mEmc.c1 * (CurrentFirstDerivative) - PreviousFirstDerivative;
      
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
      rSerializer.save("EmcParameters", mEmc);
      // rSerializer.save("StepVariable", mpStepVariable);
    };

    virtual void load(Serializer& rSerializer) override
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

#endif // KRATOS_EMC_STEP_METHOD defined
