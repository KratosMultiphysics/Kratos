//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MISES_HUBER_YIELD_CRITERION_H_INCLUDED )
#define  KRATOS_MISES_HUBER_YIELD_CRITERION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_criteria/yield_criterion.hpp"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) MisesHuberYieldCriterion : public YieldCriterion<THardeningLaw>
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;


    typedef YieldCriterion<THardeningLaw>                                BaseType;
    typedef typename BaseType::Pointer                            BaseTypePointer;
    typedef typename BaseType::PlasticDataType                    PlasticDataType;
    
    /// Pointer definition of MisesHuberYieldCriterion
    KRATOS_CLASS_POINTER_DEFINITION( MisesHuberYieldCriterion );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MisesHuberYieldCriterion() : BaseType() {}
   
    /// Copy constructor.
    MisesHuberYieldCriterion(MisesHuberYieldCriterion const& rOther) : BaseType(rOther) {}

    /// Assignment operator.
    MisesHuberYieldCriterion& operator=(MisesHuberYieldCriterion const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }
    
    /// Clone.
    virtual BaseTypePointer Clone() const override
    {
      return ( MisesHuberYieldCriterion::Pointer(new MisesHuberYieldCriterion(*this)) );
    }
    
    /// Destructor.
    virtual ~MisesHuberYieldCriterion() {}


    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Yield Condition
     */

    virtual double& CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition) override
    {
      KRATOS_TRY

      double Hardening = 0;

      const double& rStressNorm = rVariables.GetStressNorm();

      Hardening = this->mHardeningLaw.CalculateHardening(rVariables,Hardening);
<<<<<<< HEAD
		
=======
     
>>>>>>> b753b27... thermo mechanical modelling is working again
      rYieldCondition = rStressNorm - sqrt(2.0/3.0) * Hardening;
		
      return rYieldCondition;

      KRATOS_CATCH(" ")
    }
    
    /**
     * Calculate State Function
     */

    virtual double& CalculateStateFunction(const PlasticDataType& rVariables, double & rStateFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();
    
      const double& rStressNorm = rVariables.GetStressNorm();
    
      const double& rDeltaGamma = rVariables.GetDeltaInternalVariables()[0];
    
      double Hardening = 0;
		
      Hardening = this->mHardeningLaw.CalculateHardening( rVariables, Hardening );

      rStateFunction = rStressNorm - 2.0 * rMaterial.GetLameMuBar() * rDeltaGamma - sqrt(2.0/3.0) * ( Hardening );
		
      return rStateFunction;

      KRATOS_CATCH(" ")
    }
    
    /**
     * Calculate State Function derivative
     */

    virtual double& CalculateDeltaStateFunction(const PlasticDataType& rVariables, double & rDeltaStateFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      double DeltaHardening = 0;

      DeltaHardening = this->mHardeningLaw.CalculateDeltaHardening( rVariables, DeltaHardening );

      rDeltaStateFunction = 2.0 * rMaterial.GetLameMuBar() + (2.0/3.0) * DeltaHardening;
		
      return rDeltaStateFunction;

      KRATOS_CATCH(" ")
    }
    
    /**
     * Calculate Plastic Dissipation
     */

    virtual double& CalculatePlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation) override
    {
      KRATOS_TRY

      rPlasticDissipation = 0;
      return rPlasticDissipation;

      KRATOS_CATCH(" ")
    }
    
    /**
     * Calculate Plastic Dissipation derivative
     */
    
    virtual double& CalculateDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation) override
    {
      KRATOS_TRY

      rDeltaPlasticDissipation = 0;
      return rDeltaPlasticDissipation;

      KRATOS_CATCH(" ")
    }
    /**
     * Calculate Implex Plastic Dissipation
     */

    virtual double& CalculateImplexPlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation) override
    {
      KRATOS_TRY

      rPlasticDissipation = 0;
      return rPlasticDissipation;
    
      KRATOS_CATCH(" ")
    }      
    
    /**
     * Calculate Implex Plastic Dissipation derivative
     */
    
    virtual double& CalculateImplexDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation) override
    {
      KRATOS_TRY
      
      rDeltaPlasticDissipation = 0;
      return rDeltaPlasticDissipation;

      KRATOS_CATCH(" ")
    }
          
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
    virtual std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "YieldCriterion" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "MisesHuberYieldCriterion";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "MisesHuberYieldCriterion Data";
    }


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

  }; // Class MisesHuberYieldCriterion

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MISES_HUBER_YIELD_CRITERION_H_INCLUDED  defined 


