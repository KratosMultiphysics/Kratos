// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solid_mechanics_application.h"
#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.hpp"
#include "pfem_solid_mechanics_application.h"

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
   double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

   double PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);

   rStateFunction = pow(DeviatoricQ/ShearM, 2.0);
//   rStateFunction += pow(1/Beta, 2.0)*(MeanStress*(MeanStress-PreconsolidationStress) - pow(PreconsolidationStress, 2.0)*( 1.0-pow(Beta,2.0)));
   rStateFunction += pow(1.0/Beta, 2.0)*(MeanStress*(MeanStress-2.0*PreconsolidationStress)); // - pow(PreconsolidationStress, 2.0)*(1.0-pow(Beta, 2.0)));


//if (MeanStress > -10.0)
//    rStateFunction = -1.0;

//   std::cout << "P " << MeanStress << " Q " << DeviatoricQ << " PC " << PreconsolidationStress << " st " << rStateFunction << std::endl;
//   std::cout << " SV " << rStressVector << std::endl;


   rStateFunction /= 1000.0;
//   rStateFunction = -1.0;


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

/* CRITICAL STATE LINE DEPENDENCY OF LODE ANGLE
    Matrix StressMatrix;
    double DeviatoricDeterminant; //Buscar donde se calcula
    StressMatrix = MathUtils<double>::StressVectorToTensor(rStressVector);
    for (unsigned int i = 0; i<3; ++i)
       StressMatrix(i,i) -= rMeanPressure;

    DeviatoricDeterminant = 2.0;

    double LodeAngleSinus;
    LodeAngleSinus = DeviatoricDeterminant/2.0 /  pow(  rDeviatoricQ, 3.0);


    double M = 0;
   
    //CONSTANTES 
    double Mmax = 1;
    double AlphaFourth = 0.3;
  
    M = Mmax * pow( 2*AlphaFourth/  ( 1 + AlphaFourth + (1-AlphaFourth)*LodeAngleSinus), 1.0/4.0);
*/

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
    double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
  
    rYieldFunctionD = 2.0/pow(Beta, 2.0)*(MeanStress-PreconsolidationStress)*IdentityVector;
   
    rYieldFunctionD += 3.0 * ShearVector/pow(ShearM, 2.0) ;
 
    rYieldFunctionD /= 1000.0;

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
