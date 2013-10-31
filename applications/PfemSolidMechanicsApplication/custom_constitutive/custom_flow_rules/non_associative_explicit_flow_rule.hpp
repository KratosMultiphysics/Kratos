////
// NON ASSOCIATIVE, WITH HARDENING, PLASTIC FLOW RULE COMPUTED IN EXPLICITEDLY (maybe with runge-kutta)



#if !defined(KRATOS_NON_ASSOCIATIVE_EXPLICIT_PLASTIC_FLOW_RULE_H_INCLUDED)
#define      KRATOS_NON_ASSOCIATIVE_EXPLICIT_PLASTIC_FLOW_RULE_H_INCLUDED


// System includes

// External includes
#include<cmath>
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
  class NonAssociativeExplicitPlasticFlowRule
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
      KRATOS_CLASS_POINTER_DEFINITION(NonAssociativeExplicitPlasticFlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NonAssociativeExplicitPlasticFlowRule();

    /// Copy constructor.
    NonAssociativeExplicitPlasticFlowRule(const NonAssociativeExplicitPlasticFlowRule & rOther);

    /// Assignment operator.
    NonAssociativeExplicitPlasticFlowRule& operator=(NonAssociativeExplicitPlasticFlowRule const& rOther) {FlowRule::operator=(rOther); return *this;} ;
	
    /// Destructor.
    virtual ~NonAssociativeExplicitPlasticFlowRule();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    //bool CalculateReturnMapping(  RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix, const Vector& rPrincipalSrainTrial );

    //virtual bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix);

    virtual bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix &rDeltaHenckyStrain, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen);

    virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables );


    virtual void ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rLeftCauchyGreenMatrix, const double& rAlpha, Matrix& rElasticMatrix);

      virtual void ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix)
    {
	    KRATOS_ERROR(std::logic_error, "calling not the base class but another function in FlowRule ... illegal operation!!","");
    };

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

	
    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    Matrix MyCrossProduct(const Vector& rA, const Vector& rB);


    virtual void CalculateKirchhoffStressVector(const Vector& rElasticHenckyStrain, Vector& rNewStressVector)
    {
	KRATOS_ERROR(std::logic_error, "Calling the base class function in NonAss FlowRule ... illegal operation!", "");
    };

    virtual void ComputePlasticHardeningParameter(const Vector& rStressVector, const double& rAlpha, double& rH)
    {
	KRATOS_ERROR(std::logic_error, "Calling the base class function in NonAss FlowRule ... illegal operation!", "");
    };

    Matrix ConvertHenckyStrainToCauchyGreenTensor(const Vector& rElasticHenckyStrain);
   
    Vector ConvertCauchyGreenTensorToHenckyStrain(const Matrix& rCauchyGreenTensor);

    void UpdateDerivatives(const Vector& rHenckyElasticStrain, AuxiliarDerivativesStructure & rAuxiliarDerivatives, const double& EquivalentPlasticStrian) ;

    virtual void CalculatePlasticPotentialDerivatives(const Vector& rPrincipalStress, Vector& rFirstDerivative, Matrix& rSecondDerivative)
    {
	KRATOS_ERROR(std::logic_error, "Calling the base class function in NonAss FlowRule ... illegal operation!", "");
    };


    void CalculateOneExplicitStep(const Vector& rHenckyStrainIncrement, const Vector& rPreviousElasticHenckyStrain, InternalVariables& rPlasticVariables, Vector& rNewElasticHenckyStrain, Vector& rNewStressVector, double& rNewEquivalentPlasticStrain, const bool& rElastoPlasticBool, double rStressErrorMeasure);

    void CalculateExplicitSolution( const Vector& rHenckyStrainIncrement, const Vector& rPreviousElasticHenckyStrain, InternalVariables& rPlasticVariables, Vector& rNewElasticHenckyStrain, Vector& rNewStressVector,  const bool& rElastoPlasticBool, const double& rTolerance); 


    void CalculateExplicitSolutionWithChange( const Vector& rHenckyStrainIncrement, const Vector& rPreviousElasticHenckyStrain, InternalVariables& rPlasticVariables, Vector& rNewElasticHenckyStrain, Vector& rNewStressVector,  const double& rTolerance); 


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

//    virtual void save(Serializer& rSerializer) const;

//    virtual void load(Serializer& rSerializer);

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



#endif //KRATOS-NON_ASSOCIATIVE_EXPLICIT_PLASTIC_FLOW_RULE_H_INCLUDED
