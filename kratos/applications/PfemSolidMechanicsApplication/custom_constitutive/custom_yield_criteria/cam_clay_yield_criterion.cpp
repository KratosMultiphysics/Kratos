// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solid_mechanics_application.h"
#include "custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"

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
   
   const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

   double PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);

   double ThirdInvariantEffect = EvaluateThirdInvariantEffect(rStressVector);

   rStateFunction = pow(ThirdInvariantEffect * DeviatoricQ/ShearM, 2.0);
   rStateFunction += (MeanStress * (MeanStress - PreconsolidationStress) );



//   std::cout << "P " << MeanStress << " Q " << DeviatoricQ << " PC " << PreconsolidationStress << " st " << rStateFunction << std::endl;
//   std::cout << " SV " << rStressVector << std::endl;
   //rStateFunction = -1.0;
   return rStateFunction; 
}


void CamClayYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
{
    double PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);
    const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
    double MeanStress;
    double DeviatoricQ;


    this->CalculateInvariants( rStressVector, MeanStress, DeviatoricQ);
    double ThirdInvariantEffect = EvaluateThirdInvariantEffect(rStressVector);

    Vector IdentityVector = ZeroVector(6);
    for (unsigned int i = 0; i<3; ++i)
          IdentityVector(i) = 1.0/3.0;

    Vector ShearVector = ZeroVector(6);
    for (unsigned int i = 0; i<3; ++i)
       ShearVector(i) = rStressVector(i) - MeanStress;

    for (unsigned int i = 3; i<6; ++i)
       ShearVector(i) = 2.0*rStressVector(i);
 
  
    rYieldFunctionD = ( 2.0*MeanStress - PreconsolidationStress) * IdentityVector; 
 
    rYieldFunctionD += 3.0 * ShearVector * pow(ThirdInvariantEffect/ShearM, 2.0) ;

    //CalculateAndAddThirdInvDerivative( rStressVector, rYieldFunctionD);


}

double CamClayYieldCriterion::EvaluateThirdInvariantEffect( const Vector& rStressVector)
{


   double Effect = 1.0;
   double Friction = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
   if (Friction < 0.0) {
      return 1.0;
   }
   Friction *= GetPI() / 180.0;

   Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rStressVector );

   double MeanStress = 0.0;
   for (unsigned int i = 0; i < 3; ++i)
      MeanStress += StressTensor(i,i);
   MeanStress /=3;

   double J2InvSQ = 0.0;
   for (unsigned int i = 0; i < 3; ++i)
      J2InvSQ += pow(StressTensor(i,i) - MeanStress, 2.0);

   J2InvSQ += 2.0*pow(StressTensor(0,1), 2.0);
   J2InvSQ += 2.0*pow(StressTensor(0,2), 2.0);
   J2InvSQ += 2.0*pow(StressTensor(1,2), 2.0);

   J2InvSQ = sqrt( J2InvSQ/2.0);

   for (unsigned int i = 0; i < 3; ++i)
      StressTensor(i,i) -= MeanStress;
   double LodeSin = MathUtils<double>::Det(StressTensor);

   LodeSin = 3.0*sqrt(3.0)/2.0 * LodeSin / pow( J2InvSQ, 3.0);



   double epsi = 1.0e-5;
   double LodeAngle;
   if ( fabs( LodeSin ) > 1.0-epsi) {
      LodeAngle = -30.0*GetPI() / 180.0 * LodeSin / fabs(LodeSin);
   }
   else if ( J2InvSQ < 10.0*epsi) {
      LodeAngle = 30.0*GetPI() / 180.0;
   } 
   else {
      LodeAngle = std::asin( -LodeSin) / 3.0;
   }


   //std::cout << " JUST TO CHECK: rSV " << rStressVector << " LDOE " << LodeAngle/(GetPI()/180.0) << " lodeSin " << LodeSin << " invJ2 " << J2InvSQ << std::endl;

   double LodeCut = GetSmoothingLodeAngle();

   if ( fabs(LodeAngle)  < LodeCut) {

      Effect = std::cos(LodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::sin(LodeAngle); 

   }
   else {

      double A, B;
      GetSmoothingConstants(A, B, LodeAngle);
      Effect = A + B*std::sin(3.0*LodeAngle);
      //std::cout << " LODE " << LodeAngle / (GetPI() / 180.0) <<  " lodeCUT " << LodeCut/(GetPI() / 180.0) << std::endl;

   }

   //Effect /= (3.0-std::sin(Friction)) / 2.0/sqrt(3.0);
   Effect /= ( sqrt(3)/6) * (3.0 - std::sin(Friction) );

   //std::cout << " ---InTheEffect: Lode " << LodeAngle*180.0/GetPI() << " and Effect " << Effect << std::endl;
   return Effect;

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


void CamClayYieldCriterion::CalculateAndAddThirdInvDerivative(const Vector& rStressVector, Vector& rYieldFunctionD)
{

   // LAS PROPERAS LINEAS SON UN DESASTRE PORQUE ES UN COPY PASTE DE UN TROZO QUE SALE ARRIBA

   double Effect = 1.0;
   double Friction = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
   if (Friction < 0.0) {
      return;
   }
   Friction *= GetPI() / 180.0;

   Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rStressVector );

   double MeanStress = 0.0;
   for (unsigned int i = 0; i < 3; ++i)
      MeanStress += StressTensor(i,i);
   MeanStress /=3;

   double J2InvSQ = 0.0;
   for (unsigned int i = 0; i < 3; ++i)
      J2InvSQ += pow(StressTensor(i,i) - MeanStress, 2.0);

   J2InvSQ += 2.0*pow(StressTensor(0,1), 2.0);
   J2InvSQ += 2.0*pow(StressTensor(0,2), 2.0);
   J2InvSQ += 2.0*pow(StressTensor(1,2), 2.0);

   J2InvSQ = sqrt( J2InvSQ/2.0);

   for (unsigned int i = 0; i < 3; ++i)
      StressTensor(i,i) -= MeanStress;
   double LodeSin = MathUtils<double>::Det(StressTensor);

   LodeSin = 3.0*sqrt(3.0)/2.0 * LodeSin / pow( J2InvSQ, 3.0);



   double epsi = 1.0e-5;
   double LodeAngle;
   if ( fabs( LodeSin ) > 1.0-epsi) {
      LodeAngle = -30.0*GetPI() / 180.0 * LodeSin / fabs(LodeSin);
   }
   else if ( J2InvSQ < 10.0*epsi) {
      LodeAngle = 30.0*GetPI() / 180.0;
   } 
   else {
      LodeAngle = std::asin( -LodeSin) / 3.0;
   }


   double LodeCut = GetSmoothingLodeAngle();
   double C2, C3;
   double EffectDeriv;
   const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

   if ( fabs(LodeAngle)  < LodeCut) {

      Effect = std::cos(LodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::sin(LodeAngle); 
      EffectDeriv = -std::sin(LodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::cos(LodeAngle);

      C2 = -std::tan(3.0*LodeAngle) * 6.0 * Effect*EffectDeriv * pow( J2InvSQ/ShearM, 2.0);

      C3 = -6.0*sqrt(3.0) * Effect * EffectDeriv;
      C3 /= std::cos( 3.0*LodeAngle) * 2.0*  pow( ShearM, 2.0);

   }
   else {

      double A, B;
      GetSmoothingConstants(A, B, LodeAngle);
      Effect = A + B*std::sin(3.0*LodeAngle);

      EffectDeriv = 3.0*B; // this cos is ommited because in C2 is a cos and ... * cos(3.0*LodeAngle). Idem in C3;

      C2 = -std::sin(3.0*LodeAngle) *   Effect*EffectDeriv * 6.0* pow( J2InvSQ/ShearM, 2.0);

      C3 = -6.0*sqrt(3.0) * Effect * EffectDeriv;
      C3 /=  2.0* pow( ShearM, 2.0);

   }

   double Adimm = (3.0-std::sin(Friction)) * sqrt(3.0) / 6.0;
   C2 /= pow(Adimm, 2.0);
   C3 /= pow(Adimm, 2.0);


   Vector ShearStress = rStressVector;
   for (unsigned int i = 0; i < 3; ++i)
      ShearStress(i) -= MeanStress;

   Vector C2Vector = ShearStress;
   for (unsigned int i = 3; i < 6; ++i )
      C2Vector(i) *= 2.0;

   C2Vector /= 2.0 * J2InvSQ;

   Vector C3Vector = ZeroVector(6);


   // FALTER TERMES
   C3Vector(0) = ShearStress(1)*ShearStress(2) - pow( ShearStress(4), 2.0); 
   C3Vector(1) = ShearStress(2)*ShearStress(0) - pow( ShearStress(5), 2.0); 
   C3Vector(2) = ShearStress(0)*ShearStress(1) - pow( ShearStress(3), 2.0); 

   C3Vector(3) = 2.0 * ( ShearStress(4)*ShearStress(5) - ShearStress(2)*ShearStress(3));
   C3Vector(4) = 2.0 * ( ShearStress(5)*ShearStress(3) - ShearStress(0)*ShearStress(4));
   C3Vector(5) = 2.0 * ( ShearStress(3)*ShearStress(4) - ShearStress(1)*ShearStress(5));

   for (unsigned int i = 0; i < 3; ++i)
      C3Vector(i) += pow(J2InvSQ, 2.0) / 3.0;

   Vector ThisDerivative =  C2*C2Vector + C3*C3Vector;

   rYieldFunctionD += ThisDerivative;


}

//************************* YIELD FUNCTION DERIVATIVE ******************
//**********************************************************************

void CamClayYieldCriterion::GetSmoothingConstants(double& rA, double& rB, const double& rLodeAngle)
{

   
    double SmoothingAngle = this->GetSmoothingLodeAngle();
    double FrictionAngle = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
    FrictionAngle *= GetPI() / 180.0;

    double Sign = 1.0;
    if ( rLodeAngle < 0.0)
           Sign = -1.0;

    rA = 3.0 +  std::tan(SmoothingAngle) * std::tan(3.0*SmoothingAngle) + Sign * (std::tan( 3.0*SmoothingAngle) - 3.0*std::tan(SmoothingAngle)) * std::sin( FrictionAngle) / sqrt(3.0);
    rA *= (1.0/3.0) * std::cos( SmoothingAngle );

    //rSmoothingConstants.B = -1.0/ ( 3.0*std::cos( SmoothingAngle) ) * ( Sign * std::sin(SmoothingAngle) + std::sin(FrictionAngle)*std::cos(SmoothingAngle) / sqrt(3.0));
    rB = -1.0 * ( Sign* std::sin(SmoothingAngle) + std::sin(FrictionAngle)*std::cos(SmoothingAngle) / sqrt(3.0) ) / ( 3.0*std::cos(3.0*SmoothingAngle) );



}
double CamClayYieldCriterion::GetSmoothingLodeAngle()
{
    return 29.9*GetPI()/180.0;
}


double CamClayYieldCriterion::GetPI()
{
   return 3.14159265359;
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
