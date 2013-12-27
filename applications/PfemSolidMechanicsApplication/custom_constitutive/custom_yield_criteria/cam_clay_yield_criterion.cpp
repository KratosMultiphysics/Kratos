// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solid_mechanics_application.h"
#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.hpp"

namespace Kratos
{


//*******************************CONSTRUCTOR******************************************
//************************************************************************************
CamClayYieldCriterion::CamClayYieldCriterion()
	:YieldCriterion()
{
   
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

CamClayYieldCriterion::CamClayYieldCriterion(HardeningLawPointer pHardeningLaw)
	:YieldCriterion(pHardeningLaw)
{
   
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

CamClayYieldCriterion& CamClayYieldCriterion::operator=(CamClayYieldCriterion const& rOther)
{
   YieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

CamClayYieldCriterion::CamClayYieldCriterion(CamClayYieldCriterion const& rOther)
	:YieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

CamClayYieldCriterion::~CamClayYieldCriterion()
{
}



//************************* CALCULATE YIELD FUNCTION  ******************
//**********************************************************************

double& CamClayYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha)
{
 
   double MeanStress;
   double DeviatoricQ;
   this->CalculateInvariants( rStressVector, MeanStress, DeviatoricQ);
   
   double Beta = 1.0;
   double ShearM = 1.4;

   double PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);

   rStateFunction = pow(DeviatoricQ/ShearM, 2.0);
//   rStateFunction += pow(1/Beta, 2.0)*(MeanStress*(MeanStress-PreconsolidationStress) - pow(PreconsolidationStress, 2.0)*( 1.0-pow(Beta,2.0)));
   rStateFunction += pow(1.0/Beta, 2.0)*(MeanStress*(MeanStress-2.0*PreconsolidationStress) - pow(PreconsolidationStress, 2.0)*(1.0-pow(Beta, 2.0)));
if (MeanStress > -10.0)
    rStateFunction = -1.0;


//std::cout << " rSF " << rStateFunction << " PRECONSO " << PreconsolidationStress << " meanStress " << MeanStress << " devQ " << DeviatoricQ << std::endl;

   return rStateFunction; 
}



void CamClayYieldCriterion::CalculateInvariants(const Vector& rStressVector, double& rMeanPressure, double& rDeviatoricQ)
{
    rMeanPressure = 0.0;
    for (unsigned int i = 0; i<3; ++i)
         rMeanPressure += rStressVector(i);

    rMeanPressure /= 3.0;

    rDeviatoricQ = 0.0;
   
    for (unsigned int i = 0; i<3; ++i)
      rDeviatoricQ += pow( rStressVector(i) - rMeanPressure, 2.0);

    for (unsigned int i = 3; i<6; ++i)
      rDeviatoricQ += 2.0*pow(rStressVector(i), 2.0);
   
    rDeviatoricQ = pow( 3.0/2.0*rDeviatoricQ, 1.0/2.0);


}

//************************* YIELD FUNCTION DERIVATIVE ******************
//**********************************************************************
void CamClayYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
{
    double PreconsolidationStress;
    PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);
    double MeanStress;
    double DeviatoricQ;


    this->CalculateInvariants( rStressVector, MeanStress, DeviatoricQ);

    Vector IdentityVector = ZeroVector(6);
    for (unsigned int i = 0; i<3; ++i)
          IdentityVector(i) = 1.0/3.0;

    Vector ShearVector = ZeroVector(6);
    for (unsigned int i = 0; i<3; ++i)
       ShearVector(i) = rStressVector(i) - MeanStress;

    for (unsigned int i = 3; i<6; ++i)
       ShearVector(i) = 2.0*rStressVector(i);
 
    double Beta = 1.0;
    double ShearM = 1.4;
  
    rYieldFunctionD = 2.0/pow(Beta, 2.0)*(MeanStress-PreconsolidationStress)*IdentityVector;
   
    rYieldFunctionD += 3.0/pow(ShearM, 2.0) * ShearVector;

}
void CamClayYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void CamClayYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
