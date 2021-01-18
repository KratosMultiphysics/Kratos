//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_BORJA_CAM_CLAY_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_BORJA_CAM_CLAY_PLASTIC_FLOW_RULE_H_INCLUDED


// System includes

// External includes

// Project includes
//#include "custom_constitutive/custom_flow_rules/cam_clay_explicit_plastic_flow_rule.hpp"
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
  class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) BorjaCamClayExplicitFlowRule 
	  :public NonAssociativeExplicitPlasticFlowRule
  {
  


   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION(BorjaCamClayExplicitFlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BorjaCamClayExplicitFlowRule();

    /// Initialization constructor.
    BorjaCamClayExplicitFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    BorjaCamClayExplicitFlowRule(BorjaCamClayExplicitFlowRule const& rOther);

    /// Assignment operator.
    BorjaCamClayExplicitFlowRule& operator=(BorjaCamClayExplicitFlowRule const& rOther);
	
    // CLONE
    FlowRule::Pointer Clone() const override;

    /// Destructor.
    virtual ~BorjaCamClayExplicitFlowRule();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
      virtual void EvaluateMeanStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress);

      virtual void EvaluateDeviatoricStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, Vector& rDeviatoricStress);

      void ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix) override;

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
    // std::string Info() const override;

    // /// Print information about this object.
    // void PrintInfo(std::ostream& rOStream) const override;

    // /// Print object's data.
    // void PrintData(std::ostream& rOStream) const override;


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

      void CalculateKirchhoffStressVector(const Vector& rHencyStrainVector, Vector& rKirchhoffStressVector) override;


      void EvaluateMeanStress(const Vector& rHenckyStrainVector, double& rMeanStress);


      void ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH) override;

//    bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables );


    //void UpdateConfiguration( ExponentialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );	  
    
    void CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative) override; 
           
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

#endif // KRATOS_MATSUOKA_NAKAI_PLASTIC_FLOW_RULE_H_INCLUDED  defined 
