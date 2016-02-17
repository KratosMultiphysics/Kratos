//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED )
#define  KRATOS_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
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
   */
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) LinearAssociativePlasticFlowRule
	  :public NonLinearAssociativePlasticFlowRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of LinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION( LinearAssociativePlasticFlowRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LinearAssociativePlasticFlowRule();

    /// Initialization constructor.
    LinearAssociativePlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    LinearAssociativePlasticFlowRule(LinearAssociativePlasticFlowRule const& rOther);

    /// Assignment operator.
    LinearAssociativePlasticFlowRule& operator=(LinearAssociativePlasticFlowRule const& rOther);
	
    /// Destructor.
    virtual ~LinearAssociativePlasticFlowRule();


    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this flow rule
     */
    virtual FlowRule::Pointer Clone() const;

    ///@}
    ///@name Operations
    ///@{
    
    //void CalculateScalingFactors(const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors );

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    // /// Turn back information as a string.
    // virtual std::string Info() const;

    // /// Print information about this object.
    // virtual void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const;


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

    bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters );


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

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class LinearAssociativePlasticFlowRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  // /// input stream function
  // inline std::istream& operator >> (std::istream& rIStream,
  // 				    LinearAssociativePlasticFlowRule& rThis);

  // /// output stream function
  // inline std::ostream& operator << (std::ostream& rOStream,
  // 				    const LinearAssociativePlasticFlowRule& rThis)
  // {
  //   rThis.PrintInfo(rOStream);
  //   rOStream << std::endl;
  //   rThis.PrintData(rOStream);

  //   return rOStream;
  // }
  ///@}

  ///@} addtogroup block



  ///@}
  ///@ Template Operations
  ///@{


  ///@}


}  // namespace Kratos.

#endif // KRATOS_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED  defined 


