// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_flow_rules/borja_cam_clay_explicit_plastic_flow_rule.hpp"
#include "utilities/math_utils.h"
#include "includes/ublas_interface.h"
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{


// Version of cam clay non-linear elasticity with a linear shear modulus and without volumetric-shear elastic couppling

//************ CONSTRUCTOR ***********
BorjaCamClayExplicitFlowRule::BorjaCamClayExplicitFlowRule()
   :CamClayExplicitFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

BorjaCamClayExplicitFlowRule::BorjaCamClayExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
	:CamClayExplicitFlowRule(pYieldCriterion)
{
   
}

//********* ASSIGMENT OPERATOR
BorjaCamClayExplicitFlowRule& BorjaCamClayExplicitFlowRule::operator=(BorjaCamClayExplicitFlowRule const& rOther)
{
	CamClayExplicitFlowRule::operator=(rOther);
	return *this;

}



//********** COPY CONSTRUCTOR *********
BorjaCamClayExplicitFlowRule::BorjaCamClayExplicitFlowRule(BorjaCamClayExplicitFlowRule const& rOther)
      :CamClayExplicitFlowRule(rOther)
{
}

//*******   CLONE ********
FlowRule::Pointer BorjaCamClayExplicitFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new BorjaCamClayExplicitFlowRule(*this));
  return p_clone;
}



// ********** DESTRUCTOR **************
BorjaCamClayExplicitFlowRule::~BorjaCamClayExplicitFlowRule()
{
}



void BorjaCamClayExplicitFlowRule::CalculateMeanStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress)
{

    double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;    

    double DeviatoricStrain2Norm = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
	DeviatoricStrain2Norm += pow(rDeviatoricStrainVector(i), 2.0);

    for (unsigned int i = 3; i < 6; ++i)
	DeviatoricStrain2Norm += 2.0*pow(rDeviatoricStrainVector(i)/2.0, 2.0);

    rMeanStress = -ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope) * (1.0 + 1.0*AlphaShear*DeviatoricStrain2Norm/SwellingSlope);

}




void BorjaCamClayExplicitFlowRule::CalculateDeviatoricStress(const double& rVolumetricStrain, const Vector & rDeviatoricStrainVector, Vector& rDeviatoricStress)
{

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;    
    double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
    double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];
    double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];


    rDeviatoricStress = rDeviatoricStrainVector;
    double ShearModulus = AlphaShear*ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope);
    rDeviatoricStress *= 2.0*( ShearModulus + ConstantShearModulus);
    
    for (unsigned int i = 3; i<6; ++i){
         rDeviatoricStress(i) /= 2.0;  // BECAUSE VOIGT NOTATION
   }
}


void BorjaCamClayExplicitFlowRule::ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix )
{

    Matrix FourthOrderIdentity = ZeroMatrix(6);
    for (unsigned int i = 0; i<3; ++i)
       FourthOrderIdentity(i,i) = 1.0;

    for (unsigned int i = 3; i<6; ++i)
       FourthOrderIdentity(i,i) = 0.50;
// VOIGT NOTATION AND NOT KELVIN

    Matrix IdentityCross = ZeroMatrix(6);
    for (unsigned int i = 0; i<3; ++i) {
         for (unsigned int j = 0; j<3; ++j) {
            IdentityCross(i,j) = 1.0;
         }
    }


   Vector StressVector = ZeroVector(6);
   this->CalculateKirchhoffStressVector(rElasticStrainVector, StressVector);
   

   double MeanStress = 0.0;
   double VolumetricStrain = 0.0;
   for (unsigned int i = 0; i<3; i++) {
       MeanStress += StressVector(i);
       VolumetricStrain += rElasticStrainVector(i);
   }
   MeanStress /= 3.0;


   double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
   double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
   ReferencePressure /= OCR;    

   double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
   double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

   double ConstantShearModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[INITIAL_SHEAR_MODULUS];


   rElasticMatrix  = (-1.0/SwellingSlope)*MeanStress*IdentityCross;
   rElasticMatrix += 2.0*AlphaShear*ReferencePressure*std::exp(-VolumetricStrain/SwellingSlope)*(FourthOrderIdentity - (1.0/3.0)*IdentityCross);



   StressVector = rElasticStrainVector;
   for (unsigned int i = 0; i < 3; ++i)
       StressVector(i) -= VolumetricStrain / 3.0;
       

   StressVector = StressVector * 2.0*ReferencePressure * exp(-VolumetricStrain/SwellingSlope) * AlphaShear;

   for (unsigned int i = 3; i < 6 ; ++i)
         StressVector(i) /= 2.0;


   // PARTE ASQUEROSA
   for (unsigned int i = 0; i<3; ++i) {
      for (unsigned int j = 0; j<3; ++j) {
         rElasticMatrix(i,j) -= (1.0/SwellingSlope)* (StressVector(i) ); //-MeanStress);
         rElasticMatrix(i,j) -= (1.0/SwellingSlope)* (StressVector(j) ); //-MeanStress);
       }
   }

   for (unsigned int i = 0; i<3; ++i) {
      for (unsigned int j = 3; j < 6; ++j) {
         rElasticMatrix(i,j) -= 1.0*(1.0/SwellingSlope)*(StressVector(j));///2.0;
      }
   }

   for (unsigned int i = 3; i<6; ++i) {
      for (unsigned int j = 0; j<3; ++j) {
          rElasticMatrix(i,j) -= 1.0*(1.0/SwellingSlope)*(StressVector(i));///2.0;
       }
   }

  // AND THE PART DUE TO THE INITIAL SHEAR MODULUS

  rElasticMatrix +=  2.0*ConstantShearModulus * ( FourthOrderIdentity - (1.0/3.0)*IdentityCross );

}

void BorjaCamClayExplicitFlowRule::save( Serializer& rSerializer) const 
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void BorjaCamClayExplicitFlowRule::load( Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )

}

/*void LinearCamClayExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
{

   double MeanStress;
   this->CalculateMeanStress(rHenckyStrainVector, MeanStress);

   double PreconsolidationStress;
   PreconsolidationStress = mpYieldCriterion->GetHardeningLaw().CalculateHardening(PreconsolidationStress, rAlpha);



    double ReferencePreasure = 20.0;
    double SwellingSlope = 0.0078;
    double OtherSlope = 0.085;
    double Beta = 1.0;


   rH = (MeanStress-PreconsolidationStress) ;
   rH *=  (PreconsolidationStress*(1.0-pow(Beta, 2.0)) - MeanStress) ;
   rH *= PreconsolidationStress/ ( OtherSlope - SwellingSlope);
   rH *= 4.0 / pow(Beta, 4.0);



}
*/

} //end namespace kratos
