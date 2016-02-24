//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ISOTROPIC_DAMAGE_FLOW_RULE_H_INCLUDED )
#define  KRATOS_ISOTROPIC_DAMAGE_FLOW_RULE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"
#include "custom_utilities/isotropic_damage_utilities.hpp"

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) IsotropicDamageFlowRule
	  :public FlowRule
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IsotropicDamageFlowRule
      KRATOS_CLASS_POINTER_DEFINITION( IsotropicDamageFlowRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IsotropicDamageFlowRule();

    /// Initialization constructor.
    IsotropicDamageFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    IsotropicDamageFlowRule(IsotropicDamageFlowRule const& rOther);

    /// Assignment operator.
    IsotropicDamageFlowRule& operator=(IsotropicDamageFlowRule const& rOther);
	
    /// Destructor.
    virtual ~IsotropicDamageFlowRule();


    ///@}
    ///@name Operators
    ///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this flow rule
     */
    FlowRule::Pointer Clone() const;
   
    ///@}
    ///@name Operations
    ///@{

    void InitializeMaterial(YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties); 

    bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix );

    void ComputeElastoPlasticTangentMatrix( const RadialReturnVariables& rReturnMappingVariables, const Matrix& rElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticMatrix);

    bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables );

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


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

    virtual void CalculateEquivalentStrainDerivative(Vector& rEquivalentStrainDerivative, const Vector& rStrainVector, const Matrix& rLinearElasticMatrix);

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

  }; // Class IsotropicDamageFlowRule

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  ///@}

  ///@} addtogroup block



  ///@}
  ///@ Template Operations
  ///@{


  ///@}


}  // namespace Kratos.

#endif // KRATOS_ISOTROPIC_DAMAGE_FLOW_RULE_H_INCLUDED  defined 
