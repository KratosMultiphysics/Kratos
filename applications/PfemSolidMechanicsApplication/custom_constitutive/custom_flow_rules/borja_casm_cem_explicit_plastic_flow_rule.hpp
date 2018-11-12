//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_BORJA_CASM_CEM_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_BORJA_CASM_CEM_PLASTIC_FLOW_RULE_H_INCLUDED


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
  class BorjaCasmCemExplicitFlowRule 
	  :public NonAssociativeExplicitPlasticFlowRule
  {
  


  public:
    ///@name Type Definitions
    ///@{
		
    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION(BorjaCasmCemExplicitFlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BorjaCasmCemExplicitFlowRule();

    /// Initialization constructor.
    BorjaCasmCemExplicitFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    BorjaCasmCemExplicitFlowRule(BorjaCasmCemExplicitFlowRule const& rOther);

    /// Assignment operator.
    BorjaCasmCemExplicitFlowRule& operator=(BorjaCasmCemExplicitFlowRule const& rOther);
	
    // CLONE
    virtual FlowRule::Pointer Clone() const;

    /// Destructor.
    virtual ~BorjaCasmCemExplicitFlowRule();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    virtual bool CalculateReturnMappingExpl( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen);
    
    virtual void ComputeElastoPlasticTangentMatrix( const RadialReturnVariables& rReturnMappingVariables, const Matrix& rLeftCauchyGreenMatrix, const double& rAlpha, Matrix& rElasticMatrix);
    
    virtual void ComputeElasticMatrix( const Vector& rElasticStrainVector, Matrix& rElasticMatrix);
    
    virtual void EvaluateDeviatoricStress( const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, Vector& rDeviatoricStress);
    
    virtual void EvaluateMeanStress( const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress);
    
    const PlasticVariablesType& GetPlasticVariables() { return mPlasticVariables; };
    
    virtual void InitializeMaterial( YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rMaterialProperties);
    
    virtual void InitializeMaterial( const Properties& rMaterialProperties);
    
    virtual void SetPlasticVariables( const double& rInitialPreconPressure, const double& rInitialBonding); 

    //virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables);
    
    virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables);

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
		
		PlasticVariablesType mPlasticVariables;
	
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{
    
    void CalculateExplicitSolution( const Matrix& rIncrementalDeformationGradient, const Matrix& rPreviousElasticCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, const double& rTolerance);
    
    void CalculateExplicitSolutionWithChange( const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const double& rTolerance);
		
		void CalculateKirchhoffStressVector( const Vector& rHenckyStrainVector, Vector& rKirchhoffStressVector);
		
		void CalculateKirchhoffStressVector( const Matrix& rElasticLeftCauchyGreen, Vector& rStressVector);
		
		void CalculateOneExplicitPlasticStep( const Matrix& rDeltaDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, const PlasticVariablesType& rPreviousPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, PlasticVariablesType& rNewPlasticVariables, double& rDeltaPlastic);
   
		void CalculateOneExplicitStep( const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, ExplicitStressUpdateInformation& rStressUpdateInformation);
		
		void CalculatePlasticPotentialDerivatives( const Vector& rStressVector, Vector& rFirstDerivative, Matrix & rSecondDerivative, const PlasticVariablesType& rPlasticVariables);
		void CalculatePlasticPotentialDerivativesPJ2( const Vector& rStressVector, double& rFirstDerivativeP, double& rFirstDerivativeJ2, const PlasticVariablesType& rPlasticVariables);
		//void CalculatePlasticPotentialDerivativesRowe( const Vector& rStressVector, Vector& rFirstDerivative, Matrix & rSecondDerivative, const PlasticVariablesType& rPlasticVariables);
		//void CalculatePlasticPotentialDerivativesYu( const Vector& rStressVector, Vector& rFirstDerivative, Matrix & rSecondDerivative, const PlasticVariablesType& rPlasticVariables);
		
		void ComputePlasticHardeningParameter( const Vector& rHenckyStrainVector, const PlasticVariablesType& rPlasticVariables, double& rH);

		bool& EvaluateElastoPlasticUnloadingCondition( bool& rUnloadingCondition, const Matrix& rElasticLeftCauchyGreen, const Matrix& rDeltaDeformationGradient, const PlasticVariablesType& rPlasticVariables, const double& rTolerance);

		void EvaluateMeanStress( const Vector& rHenckyStrainVector, double& rMeanStress);

		void ReturnStressToYieldSurface( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rDrift, const double& rTolerance);
		//void ReturnStressToYieldSurfaceNormal( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rDrift, const double& rTolerance);
		
		void UpdateDerivatives( const Vector& rHenckyStrain, AuxiliarDerivativesStructure& rAuxiliarDerivatives, const PlasticVariablesType& rPlasticVariables);
		
		void UpdateRadialReturnVariables( RadialReturnVariables& rReturnMappingVariables, const ExplicitStressUpdateInformation& rStressUpdateInformation);

    
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
