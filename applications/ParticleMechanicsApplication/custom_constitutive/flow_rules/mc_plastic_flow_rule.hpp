//
//   Project Name:        KratosParticleMechanicsApplication $
//   Created by:          $Author:                 IIaconeta $
//   Date:                $Date:               February 2017 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_MC_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_MC_PLASTIC_FLOW_RULE_H_INCLUDED


// System includes

// External includes

#include<cmath>
// Project includes
#include "custom_constitutive/flow_rules/MPM_flow_rule.hpp"



namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

   //struct MCStressInvariants {

       //double MeanStress;
       //double J2InvSQ;
       //double LodeAngle;

   //};

   //struct MCSmoothingConstants {

       //double A;
       //double B;

   //};
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
  class MCPlasticFlowRule
          :public MPMFlowRule
  {
  


   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION( MCPlasticFlowRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MCPlasticFlowRule();

    /// Initialization constructor.
    MCPlasticFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    MCPlasticFlowRule(MCPlasticFlowRule const& rOther);

    /// Assignment operator.
    MCPlasticFlowRule& operator=(MCPlasticFlowRule const& rOther);
	
    // CLONE
    virtual MPMFlowRule::Pointer Clone() const;

    /// Destructor.
    virtual ~MCPlasticFlowRule();

	virtual bool CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen);

    virtual bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables );
    
    virtual Matrix GetElasticLeftCauchyGreen(RadialReturnVariables& rReturnMappingVariables);
    
    
    //virtual void GetPrincipalStressAndStrain(Vector& PrincipalStresses, Vector& PrincipalStrains);
    virtual void ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& alfa, Matrix& rConsistMatrix);
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    

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
    Vector mElasticPrincipalStrain;
    Vector mPlasticPrincipalStrain;
    Vector mElasticPreviousPrincipalStrain;
    Vector mPrincipalStressTrial;
    Vector mPrincipalStressUpdated;
    int mRegion;
    bool mLargeStrainBool;
    double mEquivalentPlasticStrain;
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
	void InitializeMaterial(YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw, const Properties& rProp);


      virtual void ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH);

//    virtual bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables );

	bool CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStress, Vector& rPrincipalStrain, int& region, Vector& rPrincipalStressUpdated);
    //void UpdateConfiguration( ExponentialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );	  
    
    //void CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative); 
           
	void ComputeElasticMatrix_3X3(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix);
	
	void CalculateDepSurface(Matrix& rElasticMatrix, Vector& rFNorm, Vector& rGNorm, Matrix& rAuxDep);
	
	void CalculateDepLine(Matrix& rInvD, Vector& rFNorm, Vector& rGNorm, Matrix& rAuxDep);
	
	void CalculateElastoPlasticMatrix(const RadialReturnVariables& rReturnMappingVariables, int& rRegion, Vector& DiffPrincipalStress, Matrix& rDep);

	void ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix);

    
    void CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rInverseElasticMatrix);
    void CalculateElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rElasticMatrix);
    void CalculateModificationMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rAuxT, Matrix& rInvAuxT);
    
    //void CalculateElastoPlasticMatrix(const RadialReturnVariables& rReturnMappingVariables, int& Region, Vector& DiffPrincipalStress, Matrix& Dep);
    
    void CalculateTransformationMatrix(const Matrix& rMainDirection, Matrix& rA);
    //void CalculateSmoothingConstants( MCSmoothingConstants& rSmoothingConstants, const MCStressInvariants& rStressInvariants);

   // void CalculateStressInvariants( const Vector& rStressVector, MCStressInvariants& rStressInvariants);

    double GetSmoothingLodeAngle();

    double GetPI();
 
    double GetSmoothingHiperbolic();

	//virtual void GetPrincipalStressAndStrain(Vector& PrincipalStresses, Vector& PrincipalStrains);
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
