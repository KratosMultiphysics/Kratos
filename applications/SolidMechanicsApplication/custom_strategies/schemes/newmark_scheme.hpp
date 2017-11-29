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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewmarkScheme
  {
  protected:
    
    struct NewmarkMethod
    {
      double beta;
      double gamma;

      //system constants
      double c0;
      double c1;
      double c2;
      double c3;
      double c4;
      double c5;
      double c6;
    };
    
  public:
 
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( NewmarkScheme );

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    NewmarkScheme() {}

    /// Copy Constructor.
    NewmarkScheme(NewmarkScheme& rOther) {}

    /// Clone
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

    
    void PredictVariable(NodeType& rNode, const Variable<array_1d<double, 3> >& rVariable)
    {
     KRATOS_TRY
       
     if( rVariable == DISPLACEMENT )	
       this->PredictVariable(rNode);
     else if( rVariable == VELOCITY )
       this->PredictFirstDerivative(rNode);
     else if( rVariables == ACCELERATION )
       this->PredictSecondDerivative(rNode);
     else
       KRATOS_ERROR << "This integration scheme not applies to " << rVariable << std::endl;
	
     KRATOS_CATCH( "" )
    }

    void PredictVariable(NodeType& rNode, const VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > >& rComponent)
    {
      KRATOS_TRY
	
      if( rComponent == DISPLACEMENT_X )
	this->PredictVariable(rNode,0);
      else if( rComponent == DISPLACEMENT_Y )	
	this->PredictVariable(rNode,1);
      else if( rComponent == DISPLACEMENT_Z )	
	this->PredictVariable(rNode,2);
      else if( rComponent == VELOCITY_X)
	this->PredictFirstDerivative(rNode,0);
      else if( rComponent == VELOCITY_Y )
	this->PredictFirstDerivative(rNode,1);
      else if( rComponent == VELOCITY_Z )
	this->PredictFirstDerivative(rNode,2);
      else if( rComponent == ACCELERATION_X )
	this->PredictSecondDerivative(rNode,0);    
      else if( rComponent == ACCELERATION_Y )
	this->PredictSecondDerivative(rNode,1);
      else if( rComponent == ACCELERATION_Z )
	this->PredictSecondDerivative(rNode,2);
      else
	KRATOS_ERROR << "This integration scheme not applies to " << rComponent << std::endl;      
      	
      KRATOS_CATCH( "" )
    }
    
    void UpdateVariable(NodeType& rNode, const Variable<array_1d<double, 3> >& rVariable)
    {
      KRATOS_TRY

      if( rVariable == DISPLACEMENT )	
	this->UpdateVariable(rNode);
      else if( rVariable == VELOCITY )
	this->UpdateFirstDerivative(rNode);
      else if( rVariables == ACCELERATION )
	this->UpdateSecondDerivative(rNode);
      else
	KRATOS_ERROR << "This integration scheme not applies to " << rVariable << std::endl;
	      
      KRATOS_CATCH( "" )
    }

    void UpdateVariable(NodeType& rNode, const VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > >& rComponent)
    {
      KRATOS_TRY

      if( rComponent == DISPLACEMENT_X )
	this->UpdateVariable(rNode,0);
      else if( rComponent == DISPLACEMENT_Y )	
	this->UpdateVariable(rNode,1);
      else if( rComponent == DISPLACEMENT_Z )	
	this->UpdateVariable(rNode,2);
      else if( rComponent == VELOCITY_X)
	this->UpdateFirstDerivative(rNode,0);
      else if( rComponent == VELOCITY_Y )
	this->UpdateFirstDerivative(rNode,1);
      else if( rComponent == VELOCITY_Z )
	this->UpdateFirstDerivative(rNode,2);
      else if( rComponent == ACCELERATION_X )
	this->UpdateSecondDerivative(rNode,0);    
      else if( rComponent == ACCELERATION_Y )
	this->UpdateSecondDerivative(rNode,1);
      else if( rComponent == ACCELERATION_Z )
	this->UpdateSecondDerivative(rNode,2);
      else
	KRATOS_ERROR << "This integration scheme not applies to " << rComponent << std::endl;      
      
      KRATOS_CATCH( "" )
    }

    // predict
    
    void PredictVariable(NodeType& rNode)
    {
      KRATOS_TRY
      
      KRATOS_CATCH( "" )
    }

    void PredictFirstDerivative(NodeType& rNode)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,     0);
 	          
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(CurrentVelocity) = (mNewmark.c1 * (CurrentDisplacement-PreviousDisplacement) - mNewmark.c4 * PreviousVelocity - mNewmark.c5 * PreviousAcceleration);
      
      KRATOS_CATCH( "" )      
    }

    void PredictSecondDerivative(NodeType& rNode)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION, 0);
      
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(CurrentAcceleration) = (mNewmark.c0 * (CurrentDisplacement-PreviousDisplacement) - mNewmark.c2 * PreviousVelocity - mNewmark.c3 * PreviousAcceleration);
      
      KRATOS_CATCH( "" )              
    }


    void PredictVariable(NodeType& rNode, const unsigned int c)
    {
      KRATOS_TRY
      
      KRATOS_CATCH( "" )
    }

    
    void PredictFirstDerivative(NodeType& rNode, const unsigned int c)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,     0); 
      
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      CurrentVelocity[c] =  (mNewmark.c1 * (CurrentDisplacement[c]-PreviousDisplacement[c]) - mNewmark.c4 * PreviousVelocity[c] - mNewmark.c5 * PreviousAcceleration[c]);
      
      KRATOS_CATCH( "" )      
    }

    void PredictSecondDerivative(NodeType& rNode, const unsigned int c)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION, 0);
      
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      CurrentAcceleration[c] = (mNewmark.c0 * (CurrentDisplacement[c]-PreviousDisplacement[c]) - mNewmark.c2 * PreviousVelocity[c] - mNewmark.c3 * PreviousAcceleration[c]);
      
      KRATOS_CATCH( "" )              
    }
    

    // update
    
    void UpdateVariable(NodeType& rNode)
    {
      KRATOS_TRY
      
      KRATOS_CATCH( "" )
    }

    void UpdateFirstDerivative(NodeType& rNode)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,     0);
      
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(CurrentVelocity) =  (mNewmark.c1 * (CurrentDisplacement-PreviousDisplacement) - mNewmark.c4 * PreviousVelocity - mNewmark.c5 * PreviousAcceleration);
      
      KRATOS_CATCH( "" )      
    }

    void UpdateSecondDerivative(NodeType& rNode)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION, 0);
      
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(CurrentAcceleration) = (mNewmark.c0 * (CurrentDisplacement-PreviousDisplacement) - mNewmark.c2 * PreviousVelocity - mNewmark.c3 * PreviousAcceleration);
      
      KRATOS_CATCH( "" )              
    }


    void UpdateVariable(NodeType& rNode, const unsigned int c)
    {
      KRATOS_TRY
      
      KRATOS_CATCH( "" )
    }

    
    void UpdateFirstDerivative(NodeType& rNode, const unsigned int c)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,     0);
       
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      CurrentVelocity[c] =  (mNewmark.c1 * (CurrentDisplacement[c]-PreviousDisplacement[c]) - mNewmark.c4 * PreviousVelocity[c] - mNewmark.c5 * PreviousAcceleration[c]);
      
      KRATOS_CATCH( "" )      
    }

    void UpdateSecondDerivative(NodeType& rNode, const unsigned int c)
    {
      KRATOS_TRY

      array_1d<double,3> & CurrentDisplacement  = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double,3> & CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION, 0);
      
      array_1d<double,3> & PreviousDisplacement = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double,3> & PreviousVelocity     = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double,3> & PreviousAcceleration = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      CurrentAcceleration[c] = (mNewmark.c0 * (CurrentDisplacement[c]-PreviousDisplacement[c]) - mNewmark.c2 * PreviousVelocity[c] - mNewmark.c3 * PreviousAcceleration[c]);
      
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
