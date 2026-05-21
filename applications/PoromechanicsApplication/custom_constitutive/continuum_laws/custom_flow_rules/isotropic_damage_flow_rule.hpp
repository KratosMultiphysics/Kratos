//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ISOTROPIC_DAMAGE_FLOW_RULE_H_INCLUDED)
#define  KRATOS_ISOTROPIC_DAMAGE_FLOW_RULE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/continuum_laws/custom_flow_rules/flow_rule.hpp"
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
  class KRATOS_API(POROMECHANICS_APPLICATION) IsotropicDamageFlowRule
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
    ~IsotropicDamageFlowRule() override;


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

    void InitializeMaterial(YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties) override;

    bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen) override;

    bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix ) override;

    void ComputeElastoPlasticTangentMatrix( const RadialReturnVariables& rReturnMappingVariables, const Matrix& rElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticMatrix) override;

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

    virtual bool CalculateInternalVariables(RadialReturnVariables& rReturnMappingVariables);

    virtual void CalculateEquivalentStrainDerivative(Vector& rEquivalentStrainDerivative, const RadialReturnVariables& ReturnMappingVariables, const Matrix& LinearElasticMatrix);

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
