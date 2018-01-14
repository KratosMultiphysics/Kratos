//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_NEWMARK_METHOD )
#define  KRATOS_NEWMARK_METHOD

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewmarkMethod : public TimeIntegrationMethod<TVariableType,TValueType>
  {
  protected:
    
    struct NewmarkParameters
    {
      double beta;
      double gamma;

      double delta_time;
      
      //system constants
      double c0;
      double c1;
      double c2;
      double c3;
      double c4;
      double c5;

      void SetParameters(const double& rbeta,
			 const double& rgamma,
			 const double& rdelta_time)
      {
	beta  = rbeta;
	gamma = rgamma;

	delta_time = rdelta_time;
	
	c0 = ( 1.0 / (beta * delta_time * delta_time) );
        c1 = ( gamma / (beta * delta_time) );
        c2 = ( 1.0 / (beta * delta_time) );
        c3 = ( 0.5 / (beta) - 1.0 );
        c4 = ( (gamma / beta) - 1.0  );
        c5 = ( delta_time * 0.5 * ( ( gamma / beta ) - 2.0 ) );
      }


    private:
      
      friend class Serializer;

      void save(Serializer& rSerializer) const
      {
	rSerializer.save("beta", beta);
	rSerializer.save("gamma", gamma);
	rSerializer.save("delta_time", delta_time);
	rSerializer.save("c0", c0);
	rSerializer.save("c1", c1);
	rSerializer.save("c2", c2);
	rSerializer.save("c3", c3);
	rSerializer.save("c4", c4);
	rSerializer.save("c5", c5);						
      };

      void load(Serializer& rSerializer)
      {
	rSerializer.load("beta", beta);
	rSerializer.load("gamma", gamma);
	rSerializer.load("delta_time", delta_time);
	rSerializer.load("c0", c0);
	rSerializer.load("c1", c1);
	rSerializer.load("c2", c2);
	rSerializer.load("c3", c3);
	rSerializer.load("c4", c4);
	rSerializer.load("c5", c5);
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
    
    KRATOS_CLASS_POINTER_DEFINITION( NewmarkMethod );

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    NewmarkMethod() : BaseType() {}

    /// Copy Constructor.
    NewmarkMethod(NewmarkMethod& rOther)
      :BaseType(rOther)
      ,mNewmark(rOther.mNewmark)
    {
    }

    /// Clone.
    BaseTypePointer Clone()
    {
      return BaseTypePointer( new NewmarkMethod(*this) );
    }

    /// Destructor.
    ~NewmarkMethod(){}

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
     
     double beta = 0.25;
     if (rCurrentProcessInfo.Has(NEWMARK_BETA))
       {
	 beta = rCurrentProcessInfo[NEWMARK_BETA];
       }
     double gamma = 0.5;
     if (rCurrentProcessInfo.Has(NEWMARK_GAMMA))
       {
	 gamma = rCurrentProcessInfo[NEWMARK_GAMMA];
       }
     
     mNewmark.SetParameters(beta,gamma,delta_time);
     
     KRATOS_CATCH( "" )
    }     

    // set parameters to process info
    virtual void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY
       
     rCurrentProcessInfo[NEWMARK_BETA]  = mNewmark.beta;      
     rCurrentProcessInfo[NEWMARK_GAMMA] = mNewmark.gamma;
           
     KRATOS_CATCH( "" )
    }     

    
    // get parameters
    virtual double& GetFirstDerivativeParameter(double& rParameter) override
    {
      rParameter = mNewmark.c1;
      return rParameter;
    }

    virtual double& GetSecondDerivativeParameter(double& rParameter) override
    {
      rParameter = mNewmark.c0;
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
       
     this->UpdateVariable(rNode);
     this->UpdateFirstDerivative(rNode);
     this->UpdateSecondDerivative(rNode);

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
        buffer << "NewmarkMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NewmarkMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NewmarkMethod Data";     
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

    NewmarkParameters   mNewmark;
    
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
      CurrentVariable -= this->mNewmark.delta_time * PreviousFirstDerivative;
      
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

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
       
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      CurrentVariable = PreviousVariable + this->mNewmark.delta_time * PreviousFirstDerivative + this->mNewmark.delta_time * this->mNewmark.delta_time * ( 0.5 * ( 1.0 - 2.0 * this->mNewmark.beta ) * PreviousSecondDerivative + this->mNewmark.beta * CurrentSecondDerivative );

      // variable prediction :: uniform accelerated movement **
      CurrentVariable -= this->mNewmark.delta_time * PreviousFirstDerivative + 0.5 * this->mNewmark.delta_time * this->mNewmark.delta_time * PreviousSecondDerivative;
 
      CurrentSecondDerivative = PreviousSecondDerivative;
       
      KRATOS_CATCH( "" )      
    }


    virtual void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable               = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& CurrentFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
      
      // variable prediction :: uniform accelerated movement **
      CurrentVariable += this->mNewmark.delta_time * CurrentFirstDerivative + 0.5 * this->mNewmark.delta_time * this->mNewmark.delta_time * CurrentSecondDerivative;
      
      KRATOS_CATCH( "" )
    }

    virtual void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      // consistent newmark
      CurrentFirstDerivative = (this->mNewmark.gamma/(this->mNewmark.beta*this->mNewmark.delta_time)) * (CurrentVariable-PreviousVariable) - ( (this->mNewmark.gamma/this->mNewmark.beta) - 1.0 ) * CurrentFirstDerivative - this->mNewmark.delta_time * 0.5 * ( (this->mNewmark.gamma/this->mNewmark.beta) - 2.0 ) * CurrentSecondDerivative;

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
  
  }; // Class NewmarkMethod
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{
  
  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, NewmarkMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const NewmarkMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }
  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_NEWMARK_METHOD defined
