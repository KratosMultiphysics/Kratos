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
    
    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

    typedef YieldCriterion<THardeningLaw>                                BaseType;
    typedef typename BaseType::Pointer                            BaseTypePointer;
    typedef typename BaseType::HardeningLawPointer            HardeningLawPointer;
    typedef typename BaseType::PlasticDataType                    PlasticDataType;
    
    /// Pointer definition of SimoJuYieldCriterion
    KRATOS_CLASS_POINTER_DEFINITION( SimoJuYieldCriterion );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimoJuYieldCriterion();
    
    /// Constructor.
    SimoJuYieldCriterion(HardeningLawPointer pHardeningLaw);

    /// Copy constructor.
    SimoJuYieldCriterion(SimoJuYieldCriterion const& rOther);

    /// Assignment operator.
    SimoJuYieldCriterion& operator=(SimoJuYieldCriterion const& rOther);

    /// Clone.
    virtual BaseTypePointer Clone() const override;     

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
    virtual std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "YieldCriterion" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "SimoJuYieldCriterion";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SimoJuYieldCriterion Data";
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

  }; // Class SimoJuYieldCriterion

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block


  //****************************DEFAULT CONSTRUCTOR*************************************
  //************************************************************************************
  template<class THardeningLaw>
  SimoJuYieldCriterion<THardeningLaw>::SimoJuYieldCriterion()
    :BaseType()
  {
   
  }

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************
  template<class THardeningLaw>
  SimoJuYieldCriterion<THardeningLaw>::SimoJuYieldCriterion(HardeningLawPointer pHardeningLaw)
    :BaseType(pHardeningLaw)
  {
   
  }
  
  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  template<class THardeningLaw>
  SimoJuYieldCriterion<THardeningLaw>& SimoJuYieldCriterion<THardeningLaw>::operator=(SimoJuYieldCriterion const& rOther)
  {
    BaseType::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************
  template<class THardeningLaw>
  SimoJuYieldCriterion<THardeningLaw>::SimoJuYieldCriterion(SimoJuYieldCriterion const& rOther)
    :BaseType(rOther)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  template<class THardeningLaw>
  typename YieldCriterion<THardeningLaw>::Pointer SimoJuYieldCriterion<THardeningLaw>::Clone() const
  {
    return ( SimoJuYieldCriterion::Pointer(new SimoJuYieldCriterion(*this)) );
  }

  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  template<class THardeningLaw>
  SimoJuYieldCriterion<THardeningLaw>::~SimoJuYieldCriterion()
  {
  }

  /// Operations.


  //************************** CALCULATE EQUIVALENT STRAIN *****************************
  //************************************************************************************

  template<class THardeningLaw>
  double& SimoJuYieldCriterion<THardeningLaw>::CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition)
  {
    KRATOS_TRY

    const ModelDataType& rModelData = rVariables->GetModelData();
    
    // Compute Theta parameter
    double Theta;
    
    const Matrix& rStressMatrix = rModelData.GetStressMatrix();
    
    VectorType PrincipalStresses;
    noalias(PrincipalStresses) = ConstitutiveLawUtilities::EigenValuesDirectMethod(rStressMatrix);


    double Macaulay_PrincipalStress = 0.0, Absolute_PrincipalStress = 0.0;
    
    for(unsigned int i=0; i<3; i++)
      { 
        if(PrincipalStresses[i] > 0.0)
	  {
            Macaulay_PrincipalStress += PrincipalStresses[i];
            Absolute_PrincipalStress += PrincipalStresses[i];
	  }
        else
	  {
            Absolute_PrincipalStress -= PrincipalStresses[i];
	  }
      }

    if(Absolute_PrincipalStress > 1.0e-20)
      {
        Theta = Macaulay_PrincipalStress/Absolute_PrincipalStress;
      }
    else
      {
        Theta = 0.5;
      }
    
    // Compute Equivalent Strain (rYieldCondition)
    const Matrix& rStrainMatrix = rModelData.GetStrainMatrix();
    MatrixType Auxiliar;
    noalias(Auxiliar) = prod(rStrainMatrix,rStressMatrix);
    
    double StressNorm = 0.0;
    
    for(unsigned int i=0; i<3; i++) 
      {
        StressNorm += Auxiliar(i,i);
      }
    
    const double& StrengthRatio = rModelData.GetMaterialProperties()[STRENGTH_RATIO];
    
    rYieldCondition = (Theta+(1.0-Theta)/StrengthRatio)*sqrt(StressNorm);
    
    return rYieldCondition;

    KRATOS_CATCH(" ")    
  }


  //***************************CALCULATE DAMAGE PARAMETER ******************************
  //************************************************************************************

  template<class THardeningLaw>
  double& SimoJuYieldCriterion<THardeningLaw>::CalculateStateFunction(const PlasticDataType& rVariables, double& rStateFunction)
  {
    KRATOS_TRY
        
    rStateFunction = this->mpHardeningLaw->CalculateHardening(rVariables,rStateFunction);
    
    return rStateFunction;
    
    KRATOS_CATCH(" ")    
  }


  //***************************CALCULATE DAMAGE DERIVATIVE *****************************
  //************************************************************************************

  template<class THardeningLaw>
  double& SimoJuYieldCriterion<THardeningLaw>::CalculateDeltaStateFunction(const PlasticDataType& rVariables, double& rDeltaStateFunction)
  {
    KRATOS_TRY
    
    rDeltaStateFunction = this->mpHardeningLaw->CalculateDeltaHardening(rVariables,rDeltaStateFunction);
    
    return rDeltaStateFunction;

    KRATOS_CATCH(" ")        
  }

  
}  // namespace Kratos.

#endif // KRATOS_SIMO_JU_YIELD_CRITERION_H_INCLUDED  defined 
