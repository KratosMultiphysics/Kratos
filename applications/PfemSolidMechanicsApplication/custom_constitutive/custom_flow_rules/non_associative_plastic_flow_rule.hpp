//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NON_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_NON_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED


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
  class NonAssociativePlasticFlowRule
	  :public FlowRule
  {
  public:
    ///@name Type Definitions
    ///@{
    struct AuxiliarDerivativesStructure
    { 
        Vector PlasticPotentialD;
        Vector YieldFunctionD;
        
        Matrix PlasticPotentialDD;
     };

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION( NonAssociativePlasticFlowRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonAssociativePlasticFlowRule();

    /// Initialization constructor.
    NonAssociativePlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    NonAssociativePlasticFlowRule(const NonAssociativePlasticFlowRule & rOther);

    /// Assignment operator.
    NonAssociativePlasticFlowRule& operator=(NonAssociativePlasticFlowRule const& rOther) {FlowRule::operator=(rOther); return *this;} ;
	
    /// Destructor.
    virtual ~NonAssociativePlasticFlowRule();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen);

    virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables );

    Vector GetStressVectorFromMatrix(const Matrix & rStressMatrix);

//    const Vector GetEigenValues () {return mTrialEigenValues; };
//    const Matrix GetEigenVectors() {return mEigenVectors;};
    
    void CalculatePrincipalAxisAlgorithmicTangent(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rStressMatrix, Matrix& rConsistMatrix);
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
    //Vector mTrialEigenValues;
    //Matrix mEigenVectors;
    Vector mElasticPrincipalStrain;

    bool mLargeStrainBool;
	
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeMaterial(YieldCriterionPointer& pYieldCriterion, HardeningLawPointer& pHardeningLaw, const Properties& rProp);


    void UpdateDerivatives(const Vector& rPrincipalStress, AuxiliarDerivativesStructure & rAuxiliarDerivatives);

    virtual void CalculatePlasticPotentialDerivatives(const Vector& rPrincipalStress, Vector& rFirstDerivative, Matrix& rSecondDerivative)
    {
	KRATOS_THROW_ERROR( std::logic_error, "Calling the base class function in NonAss FlowRule ... illegal operation!", "" )
    };

    void ComputePrincipalAxisStrain(RadialReturnVariables& rReturnMappingVariables, const Matrix& rStrainMatrix, Vector& rPrincipalStrain, Matrix& rEigenVectors);
//    double& CalculateNormStress ( Matrix & rStressMatrix, double& rNormStress );

	  
    bool CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, const Matrix& rInverseElasticMatrix, Vector& rPrincipalStress, Vector& rPrincipalStrain );
    
    void ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix);

    virtual Matrix GetElasticLeftCauchyGreen(const RadialReturnVariables& rReturnMappingVariables);

    void CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rInverseElasticMatrix);

    void UpdateDerivatives(const Matrix& rEigenVectors, const Matrix& rStressMatrix, AuxiliarDerivativesStructure & rAuxiliarDerivatives);

    void CalculateElasticMatrix( const RadialReturnVariables& rReturnMappingVariables, Matrix& ElasticMatrix);


//    void UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );	  

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

#endif // KRATOS_NON_ASSOCIATIVE_PLASTIC_FLOW_RULE_H_INCLUDED  defined 
