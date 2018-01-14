//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_EMC_STEP_ROTATION_METHOD )
#define  KRATOS_EMC_STEP_ROTATION_METHOD

// System includes

// External includes

// Project includes
#include "custom_strategies/time_integration_methods/emc_step_method.hpp"

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) EmcStepRotationMethod : public EmcStepMethod<TVariableType,TValueType>
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
    typedef EmcStepMethod<TVariableType,TValueType>       DerivedType;


    KRATOS_CLASS_POINTER_DEFINITION( EmcStepRotationMethod );

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    EmcStepRotationMethod() : DerivedType() {}

    /// Copy Constructor.
    EmcStepRotationMethod(EmcStepRotationMethod& rOther) : DerivedType(rOther) {}

    /// Clone.
    BaseTypePointer Clone()
    {
      return BaseTypePointer( new EmcStepRotationMethod(*this) );
    }

    /// Destructor.
    ~EmcStepRotationMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // update
    virtual void Update(NodeType& rNode) override;
     
    virtual void UpdateFirstDerivative(NodeType& rNode) override;

    virtual void UpdateSecondDerivative(NodeType& rNode) override;
       
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
        buffer << "EmcStepRotationMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "EmcStepRotationMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "EmcStepRotationMethod Data";     
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

    virtual void PredictStepVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict step variable from previous and current values
      TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);
	
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      CurrentStepVariable = CurrentVariable-PreviousVariable;
	
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
  
  }; // Class EmcStepRotationMethod
  
  ///@}

  ///@name Type Definitions
  ///@{  

  template<>
  void EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode);
 
  
  template<class TVariableType, class TValueType>
  void EmcStepRotationMethod<TVariableType,TValueType>::Update(NodeType& rNode)
  {
      KRATOS_TRY

      KRATOS_ERROR << " Calling a non compatible type update for ROTATIONS " <<std::endl;
	
      KRATOS_CATCH( "" )
  }


  template<>
  void EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::UpdateFirstDerivative(NodeType& rNode);
  
  template<class TVariableType, class TValueType>
  void EmcStepRotationMethod<TVariableType,TValueType>::UpdateFirstDerivative(NodeType& rNode)
  {
      KRATOS_TRY
	
      const TValueType& CurrentVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative        = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
 	          
      const TValueType& PreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      CurrentFirstDerivative = this->mEmc.c0 * (CurrentVariable-PreviousVariable);
      
      KRATOS_CATCH( "" )
  }

  template<>
  void EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::UpdateSecondDerivative(NodeType& rNode);
  
  template<class TVariableType, class TValueType>
  void EmcStepRotationMethod<TVariableType,TValueType>::UpdateSecondDerivative(NodeType& rNode)
  {
      KRATOS_TRY
	
      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
 	          
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      
      CurrentSecondDerivative = this->mEmc.c1 * (CurrentVariable-PreviousVariable);
	
      KRATOS_CATCH( "" )
  }
  
  ///@}
  ///@name Input and output
  ///@{
  
  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, EmcStepRotationMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const EmcStepRotationMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }


  // template<>
  // class EmcStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > : public EmcStepMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > {};
  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_EMC_STEP_ROTATION_METHOD defined
