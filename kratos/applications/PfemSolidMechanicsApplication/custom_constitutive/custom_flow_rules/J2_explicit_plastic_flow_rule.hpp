//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_J2_PLASTIC_FLOW_RULE_H_INCLUDED )
#define      KRATOS_J2_PLASTIC_FLOW_RULE_H_INCLUDED


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
  class J2ExplicitFlowRule
	  :public NonAssociativeExplicitPlasticFlowRule
  {
  


   public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NonLinearAssociativePlasticFlowRule
      KRATOS_CLASS_POINTER_DEFINITION( J2ExplicitFlowRule );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    J2ExplicitFlowRule();

    /// Initialization constructor.
    J2ExplicitFlowRule(YieldCriterionPointer pYieldCriterion);

    /// Copy constructor.
    J2ExplicitFlowRule(J2ExplicitFlowRule const& rOther);

    /// Assignment operator.
    J2ExplicitFlowRule& operator=(J2ExplicitFlowRule const& rOther);
	
    // CLONE
    virtual FlowRule::Pointer Clone() const;

    /// Destructor.
    virtual ~J2ExplicitFlowRule();


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

      virtual void CalculateKirchhoffStressVector(const Vector& rHencyStrainVector, Vector& rKirchhoffStressVector);

      virtual void ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix);

      virtual void ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH);

//    virtual bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables );


    //void UpdateConfiguration( ExponentialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );	  
    
    void CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative); 
           
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
