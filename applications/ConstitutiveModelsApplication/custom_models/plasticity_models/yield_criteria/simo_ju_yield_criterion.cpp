//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                 IPouplana  $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_criteria/simo_ju_yield_criterion.hpp"


namespace Kratos
{

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
  SimoJuYieldCriterion<THardeningLaw>::SimoJuYieldCriterion(HardeningLawType::Pointer pHardeningLaw)
    :BaseType(pHardeningLaw)
  {
   
  }
  
  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  template<class THardeningLaw>
  BaseType& SimoJuYieldCriterion<THardeningLaw>::operator=(SimoJuYieldCriterion const& rOther)
  {
    BaseType::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  SimoJuYieldCriterion<THardeningLaw>::SimoJuYieldCriterion(SimoJuYieldCriterion const& rOther)
    :BaseType(rOther)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  template<class THardeningLaw>
  BaseType::Pointer SimoJuYieldCriterion<THardeningLaw>::Clone() const
  {
    return ( BaseType::Pointer(new SimoJuYieldCriterion(*this)) );
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

    const ModelData& rModelData = rVariables->GetModelData();
    
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
    MatrixType Auxilar;
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
        
    rStateFunction = mpHardeningLaw->CalculateHardening(rVariables,rStateFunction);
    
    return rStateFunction;
    
    KRATOS_CATCH(" ")    
  }


  //***************************CALCULATE DAMAGE DERIVATIVE *****************************
  //************************************************************************************

  template<class THardeningLaw>
  double& SimoJuYieldCriterion<THardeningLaw>::CalculateDeltaStateFunction(const PlasticDataType& rVariables, double& rDeltaStateFunction)
  {
    KRATOS_TRY
    
    rDeltaStateFunction = mpHardeningLaw->CalculateDeltaHardening(rVariables,rDeltaStateFunction);
    
    return rDeltaStateFunction;

    KRATOS_CATCH(" ")        
  }



}  // namespace Kratos.
