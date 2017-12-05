//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_NEWMARK_SCHEME )
#define  KRATOS_NEWMARK_SCHEME

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "includes/variables.h"

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
  template<class TVariableType>
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewmarkScheme : public IntegrationScheme<TVariablesType>
  {
  protected:
    
    struct NewmarkMethod
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
      
    };
    
  public:
 
    ///@name Type Definitions
    ///@{

    /// BaseType
    typedef TimeIntegrationScheme<TVariableType>   BaseType;

    /// NodeType
    typedef BaseType::NodeType                     NodeType;
    
    /// KratosVariable or KratosVariableComponent    
    typedef BaseType::VariablePointer       VariablePointer;
    
    KRATOS_CLASS_POINTER_DEFINITION( NewmarkScheme );

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    NewmarkScheme() : TimeIntegrationScheme() {}

    /// Copy Constructor.
    NewmarkScheme(NewmarkScheme& rOther) : TimeIntegrationScheme(rOther) {}

    /// Clone.
    NewmarkScheme::Pointer Clone()
    {
      return NewmarkScheme::Pointer( new NewmarkScheme(*this) );
    }

    /// Destructor.
    ~NewmarkScheme(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // set parameters
    void SetParameters(const ProcessInfo& rCurrentProcessInfo)
    {
     KRATOS_TRY
       
     double delta_time = rCurrentProcessInfo[DELTA_TIME];

     if (delta_time < 1.0e-24)
        {
	  KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Scheme DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part " << std::endl;
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

    // get parameters
    double& GetFirstDerivativeParameter(double& rParameter)
    {
      rParameter = mNewmark.c1;
      return rParameter;
    }

    double& GetSecondDerivativeParameter(double& rParameter)
    {
      rParameters = mNewmark.c0
      return rParameter;
    }
    
    // predict
 
    void Predict(NodeType& rNode)
    {
     KRATOS_TRY
     
     if( mpInputVariable != nullptr ){ 

       if( *mpInputVariable == *mpVariable ){
	 UpdateFirstDerivative(rNode);
	 UpdateSecondDerivative(rNode);
       }

       if( *mpInputVariable == *mpFirstDerivative ){
	 PredictSecondDerivative(rNode);
	 PredictVariable(rNode);
       }
       
       if( *mpInputVariable == *mpSecondDerivative ){
	 PredictFirstDerivative(rNode);
	 PredictVariable(rNode);
       }
       
     }
	
     KRATOS_CATCH( "" )
    }

    void PredictVariable(NodeType& rNode)
    {
      KRATOS_TRY

      // predict variable from first and second derivative
      array_1d<double,3>& CurrentVariable          = rNode.FastGetSolutionStepValue(*mpVariable,         0);
      array_1d<double,3>& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  0);
      array_1d<double,3>& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 0);
 

      array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*mpVariable,         1);
      array_1d<double,3>& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  1);
      array_1d<double,3>& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 1);

      
      CurrentVariable = PreviousVariable + mNewmark.delta_time * PreviousFirstDerivative + mNewmark.delta_time * mNewmark.delta_time * ( 0.5 * ( 1.0 - 2.0 * mNewmark.beta ) * CurrentSecondDerivative + mNewmark.beta * CurrentSecondDerivative );
	
      KRATOS_CATCH( "" )
    }

    void PredictFirstDerivative(NodeType& rNode)
    {
      KRATOS_TRY
	
      // predict variable from second derivative
      array_1d<double,3>& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  0);
      array_1d<double,3>& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 0);
       
      array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*mpVariable,         1);
      array_1d<double,3>& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  1);
      array_1d<double,3>& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 1);

      CurrentFirstDerivative = PreviousFirstDerivative + mNewmark.delta_time * ( ( 1.0 - mNewmark.gamma ) * PreviousSecondDerivative +  mNewmark.gamma * CurrentSecondDerivative );
      
      KRATOS_CATCH( "" )      
    }

    void PredictSecondDerivative(NodeType& rNode)
    {
      KRATOS_TRY
	
      // predict second derivativ from first derivative
      array_1d<double,3>& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  0);
      array_1d<double,3>& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 0);
       
      array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*mpVariable,         1);
      array_1d<double,3>& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  1);
      array_1d<double,3>& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 1);

      CurrentSecondDerivative = ( (CurrentFirstDerivative - PreviousFirstDerivative) + mNewmark.delta_time * ( ( 1.0 - mNewmark.gamma ) * PreviousSecondDerivative ) ) * ( 1.0 / ( mNewmark.gamma * mNewmark.delta_time ) );
      
      KRATOS_CATCH( "" )      
    }
    
    // update
 
    void Update(NodeType& rNode)
    {
     KRATOS_TRY

     UpdateVariable(rNode);
     UpdateFirstDerivative(rNode);
     UpdateSecondDerivative(rNode);
	
     KRATOS_CATCH( "" )
    }
   
    void UpdateVariable(NodeType& rNode)
    {
      KRATOS_TRY
      
      KRATOS_CATCH( "" )
    }

    void UpdateFirstDerivative(NodeType& rNode)
    {
      KRATOS_TRY
	
      array_1d<double,3>& CurrentVariable          = rNode.FastGetSolutionStepValue(*mpVariable,         0);
      array_1d<double,3>& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  0);
 	          
      array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*mpVariable,         1);
      array_1d<double,3>& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  1);
      array_1d<double,3>& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 1);

      CurrentFirstDerivative = (mNewmark.c1 * (CurrentVariable-PreviousVariable) - mNewmark.c4 * PreviousFirstDerivative - mNewmark.c5 * PreviousSecondDerivative);
      
      KRATOS_CATCH( "" )      
    }

    void UpdateSecondDerivative(NodeType& rNode)
    {
      KRATOS_TRY

      array_1d<double,3>& CurrentVariable          = rNode.FastGetSolutionStepValue(*mpVariable,         0);
      array_1d<double,3>& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 0);
 	          
      array_1d<double,3>& PreviousVariable         = rNode.FastGetSolutionStepValue(*mpVariable,         1);
      array_1d<double,3>& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*mpFirstDerivative,  1);
      array_1d<double,3>& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*mpSecondDerivative, 1);
      
      CurrentSecondDerivative = (mNewmark.c0 * (CurrentVariable-PreviousVariable) - mNewmark.c2 * PreviousFirstDerivative - mNewmark.c3 * PreviousSecondDerivative);
      
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "NewmarkScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NewmarkScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "NewmarkScheme Data";     
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

    NewmarkMethod   mNewmark;
    
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
  
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@}
  
  }; // Class NewmarkScheme
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                                    NewmarkScheme& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                                    const NewmarkScheme& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream <<" : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_NEWMARK_SCHEME defined
