// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_flow_rules/cam_clay_explicit_plastic_flow_rule.hpp"
#include "utilities/math_utils.h"
#include "includes/ublas_interface.h"
namespace Kratos
{



//************ CONSTRUCTOR ***********
CamClayExplicitFlowRule::CamClayExplicitFlowRule()
   :NonAssociativeExplicitPlasticFlowRule()
{
}


//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

CamClayExplicitFlowRule::CamClayExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
	:NonAssociativeExplicitPlasticFlowRule(pYieldCriterion)
{
   
}

//********* ASSIGMENT OPERATOR
CamClayExplicitFlowRule& CamClayExplicitFlowRule::operator=(CamClayExplicitFlowRule const& rOther)
{
	NonAssociativeExplicitPlasticFlowRule::operator=(rOther);
	return *this;

}



//********** COPY CONSTRUCTOR *********
CamClayExplicitFlowRule::CamClayExplicitFlowRule(CamClayExplicitFlowRule const& rOther)
      :NonAssociativeExplicitPlasticFlowRule(rOther)
{
}

//*******   CLONE ********
FlowRule::Pointer CamClayExplicitFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new CamClayExplicitFlowRule(*this));
  return p_clone;
}



// ********** DESTRUCTOR **************
CamClayExplicitFlowRule::~CamClayExplicitFlowRule()
{
}



void CamClayExplicitFlowRule::CalculateKirchhoffStressVector(const Vector& rHenckyStrainVector, Vector& rKirchhoffStressVector)
{

    Vector IdentityVector = ZeroVector(6);

    for (unsigned int i = 0; i < 3; ++i)
       IdentityVector(i) = 1.0;

    double MeanStress; 
    double VolumetricStrain = MathUtils<double>::Dot( trans(rHenckyStrainVector), IdentityVector);

    Vector DeviatoricStrainVector; 
    DeviatoricStrainVector = rHenckyStrainVector -  (VolumetricStrain/3.0)*IdentityVector;   

    this->CalculateMeanStress(VolumetricStrain, DeviatoricStrainVector,  MeanStress );

    this->CalculateDeviatoricStress( VolumetricStrain, DeviatoricStrainVector, rKirchhoffStressVector);

    rKirchhoffStressVector += MeanStress*IdentityVector;

}

void CamClayExplicitFlowRule::CalculateMeanStress(const Vector& rHenckyStrainVector, double& rMeanStress)
{
   double VolumetricStrain = 0.0;
   for (unsigned int i = 0; i < 3; ++i)
       VolumetricStrain += rHenckyStrainVector(i);

   Vector DeviatoricStrain = rHenckyStrainVector;
  
   for (unsigned int i = 0; i < 3; ++i)
         DeviatoricStrain(i) -= VolumetricStrain/3.0;

   this->CalculateMeanStress(VolumetricStrain, DeviatoricStrain, rMeanStress);

}


void CamClayExplicitFlowRule::CalculateMeanStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress)
{

    //Reference Preasure
    //SwellingSlope
    //AlphaShear

    double ReferencePreasure = 20.0;
    double SwellingSlope = 0.0078;
    double AlphaShear = 120.0;

    double DeviatoricStrain2Norm = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
	DeviatoricStrain2Norm += pow(rDeviatoricStrainVector(i), 2.0);

    for (unsigned int i = 3; i < 6; ++i)
	DeviatoricStrain2Norm += 2.0*pow(rDeviatoricStrainVector(i)/2.0, 2.0);


    rMeanStress = -ReferencePreasure*std::exp( -rVolumetricStrain / SwellingSlope) * (1.0 + 1.0*AlphaShear*DeviatoricStrain2Norm/SwellingSlope);
 //   std::cout << " Calc Mean Stress: Volumetric " << rVolumetricStrain << " Computed P " << rMeanStress << std::endl;
}




void CamClayExplicitFlowRule::CalculateDeviatoricStress(const double& rVolumetricStrain, const Vector & rDeviatoricStrainVector, Vector& rDeviatoricStress)
{
    double ReferencePreasure = 20.0;
    double SwellingSlope = 0.0078;
    double AlphaShear = 120.0;

    rDeviatoricStress = rDeviatoricStrainVector;
    rDeviatoricStress *= 2.0*AlphaShear*ReferencePreasure* std::exp(-rVolumetricStrain / SwellingSlope);

    for (unsigned int i = 3; i<6; ++i)
         rDeviatoricStress(i) /= 2.0;  // BECAUSE VOIGT NOTATION
  
    //std::cout << "DeviatoricStrain " << rDeviatoricStrainVector << " DeviatoricStress " << rDeviatoricStress << std::endl;
    //std::cout << "ShearModulus " << 2.0*AlphaShear*ReferencePreasure*std::exp(-rVolumetricStrain/SwellingSlope) << std::endl;

}


void CamClayExplicitFlowRule::ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix )
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


   double SwellingSlope = 0.0078;
   double AlphaShear = 120.0;
   double ReferencePreasure = 20.0;


   rElasticMatrix  = (-1.0/SwellingSlope)*MeanStress*IdentityCross;
   rElasticMatrix += 2.0*AlphaShear*ReferencePreasure*std::exp(-VolumetricStrain/SwellingSlope)*(FourthOrderIdentity - (1.0/3.0)*IdentityCross);


   // PARTE ASQUEROSA
   for (unsigned int i = 0; i<3; ++i) {
      for (unsigned int j = 0; j<3; ++j) {
         rElasticMatrix(i,j) -= (1.0/SwellingSlope)* (StressVector(i)-MeanStress);
         rElasticMatrix(i,j) -= (1.0/SwellingSlope)* (StressVector(j)-MeanStress);
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

   if ( fabs(VolumetricStrain) > 0.50) {
     for (unsigned int i = 0; i<6; ++i) {
       rElasticMatrix(i,i) += 0000.1;
     }
   }
//std::cout << "ELASTIC MATRIX " << rElasticMatrix << "ElasticVector " << rElasticStrainVector << " StressV " << StressVector << "Mean Stress " << MeanStress << std::endl;

}


void CamClayExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
{

   double MeanStress;
   this->CalculateMeanStress(rHenckyStrainVector, MeanStress);

   double PreconsolidationStress;
   PreconsolidationStress = mpYieldCriterion->GetHardeningLaw().CalculateHardening(PreconsolidationStress, rAlpha);



    double ReferencePreasure = 20.0;
    double SwellingSlope = 0.0078;
    double AlphaShear = 120.0;
    double OtherSlope = 0.085;
    double Beta = 1.0;
    double VoidRatio = 0.6;


   rH = (MeanStress-PreconsolidationStress) ;
   rH *=  (PreconsolidationStress*(1.0-pow(Beta, 2.0)) - MeanStress) ;
   rH *= PreconsolidationStress/ ( OtherSlope - SwellingSlope);
   rH *= 4.0 / pow(Beta, 4.0);



}


} //end namespace kratos
