//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  IPouplana $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                December 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_criteria/modified_mises_yield_criterion.hpp"

namespace Kratos
{

  //****************************DEFAULT CONSTRUCTOR*************************************
  //************************************************************************************

  template<class THardeningLaw>
  ModifiedMisesYieldCriterion<THardeningLaw>::ModifiedMisesYieldCriterion()
    :BaseType()
  {
   
  }

  //*******************************CONSTRUCTOR******************************************
  //************************************************************************************

  template<class THardeningLaw>
  ModifiedMisesYieldCriterion<THardeningLaw>::ModifiedMisesYieldCriterion(HardeningLawType::Pointer pHardeningLaw)
    :BaseType(pHardeningLaw)
  {
   
  }
  
  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  template<class THardeningLaw>
  ModifiedMisesYieldCriterion& ModifiedMisesYieldCriterion<THardeningLaw>::operator=(ModifiedMisesYieldCriterion const& rOther)
  {
    YieldCriterion::operator=(rOther);
    return *this;
  }

  //*******************************COPY CONSTRUCTOR*************************************
  //************************************************************************************

  template<class THardeningLaw>
  ModifiedMisesYieldCriterion<THardeningLaw>::ModifiedMisesYieldCriterion(ModifiedMisesYieldCriterion const& rOther)
    :BaseType(rOther)
  {

  }

  //********************************CLONE***********************************************
  //************************************************************************************

  template<class THardeningLaw>
  BaseTye::Pointer ModifiedMisesYieldCriterion<THardeningLaw>::Clone() const
  {
    return ( BaseType::Pointer(new ModifiedMisesYieldCriterion(*this)) );
  }

  //********************************DESTRUCTOR******************************************
  //************************************************************************************

  template<class THardeningLaw>
  ModifiedMisesYieldCriterion<THardeningLaw>::~ModifiedMisesYieldCriterion()
  {
  }

  /// Operations.


  //************************** CALCULATE EQUIVALENT STRAIN *****************************
  //************************************************************************************

  template<class THardeningLaw>
  double& ModifiedMisesYieldCriterion<THardeningLaw>::CalculateYieldCondition(const PlasticDataType& rVariables, double& rYieldCondition)
  {
    KRATOS_TRY

    const ModelData& rModelData = rVariables->GetModelData();

    // Compute I1
    const MatrixType& rStrainMatrix = rModelData.GetStrainMatrix();
    double I1 = 0.0;
    
    for(unsigned int i=0; i<3; i++)
      {
        I1 += rStrainMatrix(i,i);
      }
    
    // Compute J2
    MatrixType DeviatoricStrain;
    noalias(DeviatoricStrain) = rStrainMatrix;
    
    for(unsigned int i=0; i<3; i++)
      {
        DeviatoricStrain(i,i) -= I1/3.0;
      }
    
    MatrixType Auxiliar;
    noalias(Auxiliar) = prod(DeviatoricStrain,DeviatoricStrain);
    double J2 = 0.0;
    
    for(unsigned int i=0; i<3; i++)
      {
        J2 += Auxiliar(i,i);
      }
    
    J2 *= 0.5;
    
    // Compute Equivalent Strain (rYieldCondition)
    const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
    const double& StrengthRatio = rMaterialProperties[STRENGTH_RATIO];
    const double& PoissonRatio  = rMaterialProperties[POISSON_RATIO];
        
    rYieldCondition  =  I1*(StrengthRatio-1.0)/(2.0*StrengthRatio*(1.0-2.0*PoissonRatio));
    rYieldCondition +=  sqrt( I1*I1*(StrengthRatio-1.0)*(StrengthRatio-1.0)/((1.0-2.0*PoissonRatio)*(1.0-2.0*PoissonRatio)) + J2*12.0*StrengthRatio/((1.0+PoissonRatio)*(1.0+PoissonRatio)) )/(2.0*StrengthRatio);
    
    return rYieldCondition;

    KRATOS_CATCH(" ")

  }


  //***************************CALCULATE DAMAGE PARAMETER ******************************
  //************************************************************************************

  template<class THardeningLaw>
  double& ModifiedMisesYieldCriterion<THardeningLaw>::CalculateStateFunction(const PlasticDataType& rVariables, double& rStateFunction)
  {    
    KRATOS_TRY
      
    rStateFunction = mpHardeningLaw->CalculateHardening(rVariables,rStateFunction);
    
    return rStateFunction;

    KRATOS_CATCH(" ")
  }


  //***************************CALCULATE DAMAGE DERIVATIVE *****************************
  //************************************************************************************

  template<class THardeningLaw>
  double& ModifiedMisesYieldCriterion<THardeningLaw>::CalculateDeltaStateFunction(const PlasticDataType& rVariables, double& rDeltaStateFunction)
  {
    KRATOS_TRY
          
    rDeltaStateFunction = mpHardeningLaw->CalculateDeltaHardening(rVariables,rDeltaStateFunction);
    
    return rDeltaStateFunction;

    KRATOS_CATCH(" ")
  }





}  // namespace Kratos.
