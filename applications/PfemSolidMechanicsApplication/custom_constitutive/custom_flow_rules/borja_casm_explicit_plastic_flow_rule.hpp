//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_BORJA_CASM_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_BORJA_CASM_PLASTIC_FLOW_RULE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"

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
  class BorjaCasmExplicitFlowRule 
	  :public NonAssociativeExplicitPlasticFlowRule
  {
  


  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION(BorjaCasmExplicitFlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BorjaCasmExplicitFlowRule();

    /// Initialization constructor.
    BorjaCasmExplicitFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    BorjaCasmExplicitFlowRule(BorjaCasmExplicitFlowRule const& rOther);

    /// Assignment operator.
    BorjaCasmExplicitFlowRule& operator=(BorjaCasmExplicitFlowRule const& rOther);
	
    // CLONE
    virtual FlowRule::Pointer Clone() const;

    /// Destructor.
    virtual ~BorjaCasmExplicitFlowRule();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    virtual void EvaluateMeanStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress);

    virtual void EvaluateDeviatoricStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, Vector& rDeviatoricStress);

    virtual void ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix);

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

		void CalculateKirchhoffStressVector(const Vector& rHenckyStrainVector, Vector& rKirchhoffStressVector);

		void EvaluateMeanStress(const Vector& rHenckyStrainVector, double& rMeanStress);

		virtual void ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH);

//    virtual bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables );

//		void UpdateConfiguration( ExponentialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );	  
    
		void CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative, const double& rAlpha); 

//    void CalculateInvariantsAndDerivatives(const Vector& rStressVector, InvariantsStructure& rInv);
    
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

#endif // KRATOS_MATSUOKA_NAKAI_PLASTIC_FLOW_RULE_H_INCLUDED  defined 
