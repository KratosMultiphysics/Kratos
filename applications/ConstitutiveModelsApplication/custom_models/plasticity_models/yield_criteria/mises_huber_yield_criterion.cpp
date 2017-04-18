//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_criteria/mises_huber_yield_criterion.hpp"

namespace Kratos
{

  //****************************DEFAULT CONSTRUCTOR*************************************
  //************************************************************************************

  template<class THardeningLaw>
  MisesHuberYieldCriterion<THardeningLaw>::MisesHuberYieldCriterion()
    :BaseType()
  {
   
  }

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  template<class THardeningLaw>
  MisesHuberYieldCriterion<THardeningLaw>::MisesHuberYieldCriterion(HardeningLawType::Pointer pHardeningLaw)
    :BaseType(pHardeningLaw)
  {
   
  }
  
  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  template<class THardeningLaw>
  MisesHuberYieldCriterion& MisesHuberYieldCriterion<THardeningLaw>::operator=(MisesHuberYieldCriterion const& rOther)
  {
    BaseType::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  template<class THardeningLaw>
  MisesHuberYieldCriterion<THardeningLaw>::MisesHuberYieldCriterion(MisesHuberYieldCriterion const& rOther)
    :BaseType(rOther)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  template<class THardeningLaw>
  BaseType::Pointer MisesHuberYieldCriterion::Clone() const
  {
    return ( BaseType::Pointer(new MisesHuberYieldCriterion(*this)) );
  }

  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  template<class THardeningLaw>
  MisesHuberYieldCriterion<THardeningLaw>::~MisesHuberYieldCriterion()
  {
  }

  /// Operations.


  //***************************CALCULATE YIELD CONDITION********************************
  //************************************************************************************

  template<class THardeningLaw>
  double& MisesHuberYieldCriterion<THardeningLaw>::CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition)
  {
    KRATOS_TRY

    double Hardening = 0;

    const double& rStressNorm = rVariables.GetStressNorm();

    Hardening = mpHardeningLaw->CalculateHardening(rVariables,Hardening);
		
    rYieldCondition = rStressNorm - sqrt(2.0/3.0) * Hardening;
		
    return rYieldCondition;

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE STATE FUNCTION ********************************
  //************************************************************************************

  template<class THardeningLaw>
  double& MisesHuberYieldCriterion<THardeningLaw>::CalculateStateFunction(const PlasticDataType& rVariables, double & rStateFunction)
  {
    KRATOS_TRY

    const MaterialData& rMaterial = rVariables.GetMaterialParameters();
    
    const double& rStressNorm = rVariables.GetStressNorm();
    
    const double& rDeltaGamma = rVariables.GetDeltaInternalVariables()[0];
    
    double Hardening = 0;
		
    Hardening = mpHardeningLaw->CalculateHardening( Hardening );

    rStateFunction = rStressNorm - 2.0 * rMaterial.GetLameMuBar() * rDeltaGamma - sqrt(2.0/3.0) * ( Hardening );
		
    return rStateFunction;

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE STATE FUNCTION ********************************
  //************************************************************************************

  template<class THardeningLaw>
  double& MisesHuberYieldCriterion<THardeningLaw>::CalculateDeltaStateFunction(const PlasticDataType& rVariables, double & rDeltaStateFunction)
  {
    KRATOS_TRY

    const MaterialData& rMaterial = rVariables.GetMaterialParameters();

    double DeltaHardening = 0;

    DeltaHardening = mpHardeningLaw->CalculateDeltaHardening( DeltaHardening );

    //std::cout<<" DeltaHardening "<<DeltaHardening<<std::endl;

    rDeltaStateFunction = 2.0 * rMaterial.GetLameMuBar() + (2.0/3.0) * DeltaHardening;
		
    return rDeltaStateFunction;

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE PLASTIC DISSIPATION****************************
  //************************************************************************************

  template<class THardeningLaw>
  double& MisesHuberYieldCriterion<THardeningLaw>::CalculatePlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation)
  {
    KRATOS_TRY

    rPlasticDissipation = 0;
    return rPlasticDissipation;

    KRATOS_CATCH(" ")
  }


  //**********************CALCULATE DELTA PLASTIC DISSIPATION***************************
  //************************************************************************************

  template<class THardeningLaw>
  double& MisesHuberYieldCriterion<THardeningLaw>::CalculateDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation)
  {
    KRATOS_TRY

    rDeltaPlasticDissipation = 0;
    return rDeltaPlasticDissipation;

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE PLASTIC DISSIPATION****************************
  //************************************************************************************

  template<class THardeningLaw>
  double& MisesHuberYieldCriterion<THardeningLaw>::CalculateImplexPlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation)
  {
    KRATOS_TRY

    rPlasticDissipation = 0;
    return rPlasticDissipation;
    
    KRATOS_CATCH(" ")
  }


  //**********************CALCULATE DELTA PLASTIC DISSIPATION***************************
  //************************************************************************************

  template<class THardeningLaw>
  double& MisesHuberYieldCriterion<THardeningLaw>::CalculateImplexDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation)
  {
    KRATOS_TRY
      
    rDeltaPlasticDissipation = 0;
    return rDeltaPlasticDissipation;

    KRATOS_CATCH(" ")
  }


}  // namespace Kratos.
