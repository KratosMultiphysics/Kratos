// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/linear_cam_clay_explicit_plastic_flow_rule.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


// Version of cam clay non-linear elasticity with a linear shear modulus and without volumetric-shear elastic couppling

//************ CONSTRUCTOR ***********
LinearCamClayExplicitFlowRule::LinearCamClayExplicitFlowRule()
   :CamClayExplicitFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

LinearCamClayExplicitFlowRule::LinearCamClayExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
	:CamClayExplicitFlowRule(pYieldCriterion)
{
   
}

//********* ASSIGMENT OPERATOR
LinearCamClayExplicitFlowRule& LinearCamClayExplicitFlowRule::operator=(LinearCamClayExplicitFlowRule const& rOther)
{
	CamClayExplicitFlowRule::operator=(rOther);
	return *this;

}



//********** COPY CONSTRUCTOR *********
LinearCamClayExplicitFlowRule::LinearCamClayExplicitFlowRule(LinearCamClayExplicitFlowRule const& rOther)
      :CamClayExplicitFlowRule(rOther)
{
}

//*******   CLONE ********
FlowRule::Pointer LinearCamClayExplicitFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new LinearCamClayExplicitFlowRule(*this));
  return p_clone;
}



// ********** DESTRUCTOR **************
LinearCamClayExplicitFlowRule::~LinearCamClayExplicitFlowRule()
{
}



/*void LinearCamClayExplicitFlowRule::CalculateKirchhoffStressVector(const Vector& rHenckyStrainVector, Vector& rKirchhoffStressVector)
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

void LinearCamClayExplicitFlowRule::CalculateMeanStress(const Vector& rHenckyStrainVector, double& rMeanStress)
{
   double VolumetricStrain = 0.0;
   for (unsigned int i = 0; i < 3; ++i)
       VolumetricStrain += rHenckyStrainVector(i);

   Vector DeviatoricStrain = rHenckyStrainVector;
  
   for (unsigned int i = 0; i < 3; ++i)
         DeviatoricStrain(i) -= VolumetricStrain/3.0;

   this->CalculateMeanStress(VolumetricStrain, DeviatoricStrain, rMeanStress);

}

*/
void LinearCamClayExplicitFlowRule::CalculateMeanStress(const double& rVolumetricStrain, const Vector& rDeviatoricStrainVector, double& rMeanStress)
{

    double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
    double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
    ReferencePressure /= OCR;    
    double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];


    rMeanStress = -ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope) ;
}




void LinearCamClayExplicitFlowRule::CalculateDeviatoricStress(const double& rVolumetricStrain, const Vector & rDeviatoricStrainVector, Vector& rDeviatoricStress)
{
    double YoungModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    double PoissonCoef  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    double ShearModulus = YoungModulus / 2.0 / (1.0 + PoissonCoef);
    //ShearModulus = 1200.0;
    ShearModulus = 23.5*50.0;

    rDeviatoricStress = rDeviatoricStrainVector;
    rDeviatoricStress *= 2.0*ShearModulus;

    for (unsigned int i = 3; i<6; ++i)
         rDeviatoricStress(i) /= 2.0;  // BECAUSE VOIGT NOTATION
  
    //std::cout << "DeviatoricStrain " << rDeviatoricStrainVector << " DeviatoricStress " << rDeviatoricStress << std::endl;
    //std::cout << "ShearModulus " << 2.0*AlphaShear*ReferencePreasure*std::exp(-rVolumetricStrain/SwellingSlope) << std::endl;


}


void LinearCamClayExplicitFlowRule::ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix )
{

    double YoungModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    double PoissonCoef  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
    double ShearModulus = YoungModulus / 2.0 / (1.0 + PoissonCoef);
    ShearModulus = 23.50*50.0;
    double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];

    Matrix FourthOrderIdentity = ZeroMatrix(6);
    for (unsigned int i = 0; i<3; ++i)
       FourthOrderIdentity(i,i) = 1.0;

    for (unsigned int i = 3; i<6; ++i)
      FourthOrderIdentity(i,i) = 0.5;

    Matrix IdentityCross = ZeroMatrix(6);
    for (unsigned int i = 0; i<3; ++i) {
         for (unsigned int j = 0; j<3; ++j) {
            IdentityCross(i,j) = 1.0;
         }
    }


   Vector StressVector = ZeroVector(6);
   this->CalculateKirchhoffStressVector(rElasticStrainVector, StressVector);
   

   double MeanStress = 0.0;
   for (unsigned int i = 0; i<3; i++) {
       MeanStress += StressVector(i);
   }
   MeanStress /= 3.0;

   rElasticMatrix  = (-1.0/SwellingSlope)*MeanStress*IdentityCross;
   rElasticMatrix += 2.0*ShearModulus*(FourthOrderIdentity - IdentityCross/3.0);


}

void LinearCamClayExplicitFlowRule::save( Serializer& rSerializer) const 
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void LinearCamClayExplicitFlowRule::load( Serializer& rSerializer)
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
