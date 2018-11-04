//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_NON_LINEAR_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED)
#define  KRATOS_NON_LINEAR_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) NonLinearAssociativePlasticFlowRule
	  :public FlowRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION( NonLinearAssociativePlasticFlowRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonLinearAssociativePlasticFlowRule();

    /// Initialization constructor.
    NonLinearAssociativePlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    NonLinearAssociativePlasticFlowRule(NonLinearAssociativePlasticFlowRule const& rOther);

    /// Assignment operator.
    NonLinearAssociativePlasticFlowRule& operator=(NonLinearAssociativePlasticFlowRule const& rOther);

    /// Destructor.
    ~NonLinearAssociativePlasticFlowRule() override;


    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this flow rule
     */
    FlowRule::Pointer Clone() const override;

    ///@}
    ///@name Operations
    ///@{

    bool CalculateReturnMapping(  RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix ) override;

    void CalculateScalingFactors( const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors ) override;

    bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables ) override;


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

    double& CalculateStressNorm ( Matrix & rStressMatrix, double& rStressNorm ) override;


    virtual void SetCriterionParameters( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters );


    virtual bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters);


    void UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );


    void CalculateThermalDissipation( YieldCriterion::Parameters& rCriterionParameters, ThermalVariables& rThermalVariables );


    //implex protected methods

    virtual void CalculateImplexReturnMapping( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters, Matrix& rIsoStressMatrix );

    void CalculateImplexThermalDissipation( YieldCriterion::Parameters& rCriterionParameters );


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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class NonLinearAssociativePlasticFlowRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  // inline std::istream& operator >> (std::istream& rIStream,
  // 				    NonLinearAssociativePlasticFlowRule& rThis);

  // /// output stream function
  // inline std::ostream& operator << (std::ostream& rOStream,
  // 				    const NonLinearAssociativePlasticFlowRule& rThis)
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

#endif // KRATOS_NON_LINEAR_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED  defined


