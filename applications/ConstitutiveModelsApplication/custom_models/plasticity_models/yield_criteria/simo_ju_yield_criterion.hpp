//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  IPouplana $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SIMO_JU_YIELD_CRITERION_H_INCLUDED )
#define  KRATOS_SIMO_JU_YIELD_CRITERION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_criteria/yield_criterion.hpp"

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
  template<class THardeningLaw>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SimoJuYieldCriterion
    : public YieldCriterion<THardeningLaw>
  {
  public:
    ///@name Type Definitions
    ///@{
    typedef YieldCriterion<THardeningLaw>     BaseType;

    /// Pointer definition of SimoJuYieldCriterion
    KRATOS_CLASS_POINTER_DEFINITION( SimoJuYieldCriterion );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimoJuYieldCriterion();
    
    /// Constructor.
    SimoJuYieldCriterion(HardeningLawType::Pointer pHardeningLaw);

    /// Copy constructor.
    SimoJuYieldCriterion(SimoJuYieldCriterion const& rOther);

    /// Assignment operator.
    SimoJuYieldCriterion& operator=(SimoJuYieldCriterion const& rOther);

    /// Clone.
    virtual BaseType::Pointer Clone() const override;     

    /// Destructor.
    virtual ~SimoJuYieldCriterion();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Yield Condition
     */

    double& CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition) override;

    /**
     * Calculate State Function
     */

    double& CalculateStateFunction(const PlasticDataType& rVariables, double & rStateFunction) override;

    /**
     * Calculate State Function derivative
     */

    double& CalculateDeltaStateFunction(const PlasticDataType& rVariables, double & rDeltaStateFunction) override;
    

    ///@}
    ///@name Access
    ///@{
        

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    
    /// Turn back information as a string.
    virtual std::string Info() const
    {
      std::stringstream buffer;
      buffer << "YieldCriterion" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SimoJuYieldCriterion";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}
    
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

    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    virtual void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class SimoJuYieldCriterion

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SIMO_JU_YIELD_CRITERION_H_INCLUDED  defined 
