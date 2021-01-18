//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/stress_invariants_utilities.hpp"
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

      double MeanStress, LodeAngle;
      double DeviatoricQ; // == sqrt(3)*J2

      StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, DeviatoricQ, LodeAngle);
      DeviatoricQ *= sqrt(3.0);


      const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];

      double PreconsolidationStress = 0.0;
      PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);
      double ThirdInvariantEffect = EvaluateThirdInvariantEffect(LodeAngle);

      rStateFunction = pow(ThirdInvariantEffect * DeviatoricQ/ShearM, 2);
      rStateFunction += (MeanStress * (MeanStress - PreconsolidationStress) );

      return rStateFunction; 
   }


   //*******************************CALCULATE YIELD FUNCTION DERIVATIVE *****************
   //************************************************************************************
   void CamClayYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
   {
      double PreconsolidationStress = 0.0;
      PreconsolidationStress = mpHardeningLaw->CalculateHardening(PreconsolidationStress, rAlpha);
      const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];
      double MeanStress, J2, LodeAngle;

      Vector V1, V2;

      StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
      StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2);

      double ThirdInvariantEffect = EvaluateThirdInvariantEffect( LodeAngle);


      rYieldFunctionD = ( 2.0*MeanStress - PreconsolidationStress) * V1 + 2.0 * 3.0 * pow( ThirdInvariantEffect / ShearM, 2) * J2 * V2;

      CalculateAndAddThirdInvDerivative( rStressVector, rYieldFunctionD);
   }

   //*******************************Evaluate Effect THird Invariant *** *****************
   //************************************************************************************
   double CamClayYieldCriterion::EvaluateThirdInvariantEffect( const double& rLodeAngle)
   {


      double Effect = 1.0;
      double Friction = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
      if (Friction < 1.0E-3) {
         return 1.0;
      }
      Friction *= GetPI() / 180.0;


      if ( true) {
         double LodeCut = GetSmoothingLodeAngle();

         if ( fabs( rLodeAngle)  < LodeCut) {
            Effect = std::cos( rLodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::sin(rLodeAngle); 
         }
         else {

            double A, B;
            GetSmoothingConstants(A, B, rLodeAngle);
            Effect = A + B*std::sin(3.0*rLodeAngle);
         }

         Effect /= ( sqrt(3)/6) * (3.0 - std::sin(Friction) );
      }
      else {
         double Betta = 0.9;
         Betta *= -1.0; // different lode angle definition, but this works.
         double GammaL = 6.0 / GetPI()  * std::atan(  std::sin( Friction) / sqrt(3.0) );
         double Gamma = ( 1- GammaL);
         double Alpha = 1.0 / std::cos( ( 1.0 + GammaL ) * GetPI() / 6.0 );

         Effect = Alpha * std::cos(   std::acos( Betta * std::sin(3.0*rLodeAngle ) ) / 3.0 - Gamma * GetPI() / 6.0 );

      }
      //std::cout << "LODE: " << rLodeAngle << " " << rLodeAngle * 180.0 / GetPI() << "  " << EffectPrev << " " << Effect << std::endl;



      return Effect;

   }


   //*******************************Add  derivative of Effect THird Invariant *** *****************
   //************************************************************************************
   void CamClayYieldCriterion::CalculateAndAddThirdInvDerivative(const Vector& rStressVector, Vector& rYieldFunctionD)
   {

      double Effect = 1.0;
      double Friction = this->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];
      if (Friction < 1.0E-3) {
         return ;
      }
      Friction *= GetPI() / 180.0;


      double MeanStress, J2, LodeAngle;
      Vector V1, V2, V3;

      StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2, LodeAngle);
      // since I will be dividing by J2
      if ( J2 < 1E-5)
         return;

      StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2, V3);
      double C2, C3;
      double EffectDeriv;
      const double ShearM = this->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE];


      if (true) {
         double LodeCut = GetSmoothingLodeAngle();

         if ( fabs(LodeAngle)  < LodeCut) {

            Effect = std::cos(LodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::sin(LodeAngle); 
            EffectDeriv = -std::sin(LodeAngle) - 1.0/sqrt(3.0) * std::sin(Friction) * std::cos(LodeAngle);

            C2 = -std::tan(3.0*LodeAngle) * 6.0 * Effect*EffectDeriv * pow( J2/ShearM, 2);
            C2 = -6.0 * J2 / pow( ShearM, 2) * Effect * EffectDeriv * std::tan(3.0*LodeAngle);

            C3 = -6.0*sqrt(3.0) * Effect * EffectDeriv;
            C3 /= std::cos( 3.0*LodeAngle) * 2.0*  pow( ShearM, 2);

            C3 = -3 * sqrt(3.0) * Effect * EffectDeriv / ( pow( ShearM, 2) * std::cos( 3.0*LodeAngle) );
            C3 /= J2;
         }
         else {

            double A, B;
            GetSmoothingConstants(A, B, LodeAngle);
            Effect = A + B*std::sin(3.0*LodeAngle);

            C2 = - 18.0 * J2 / pow(ShearM, 2) * B * Effect * sin( 3.0*LodeAngle);
            C2 = -18.0 * J2 / pow(ShearM, 2) * B * Effect * sin(3.0*LodeAngle);

            C3 = -9.0*sqrt(3.0) * Effect * B;
            C3 /=   pow( ShearM, 2);
            C3 /= J2;

         }

         double Adimm = (3.0-std::sin(Friction)) * sqrt(3.0) / 6.0;
         C2 /= pow(Adimm, 2);
         C3 /= pow(Adimm, 2);
      }
      else {

         double Betta = 0.90;
         Betta *= -1.0; // different lode angle definition, but this works.
         double GammaL = 6.0 / GetPI()  * std::atan(  std::sin( Friction) / sqrt(3.0) );
         double Gamma = ( 1- GammaL);
         double Alpha = 1.0 / std::cos( ( 1.0 + GammaL ) * GetPI() / 6.0 );

         Effect = Alpha * std::cos(   std::acos( Betta * std::sin(3.0*LodeAngle ) ) / 3.0 - Gamma * GetPI() / 6.0 );
         EffectDeriv = Alpha * Betta * std::sin(   std::acos( Betta * std::sin(3.0*LodeAngle ) ) / 3.0 - Gamma * GetPI() / 6.0 );
         EffectDeriv /=  sqrt(  1.0 - pow( Betta * std::sin( 3.0 * LodeAngle), 2 ) );
         //EffectDeriv *= std::cos( 3.0 * LodeAngle) ; // removing cos

         //C2 = -6.0 * J2 / pow( ShearM, 2.0) * Effect * EffectDeriv * std::tan(3.0*LodeAngle); // removing cos
         C2 = -6.0 * J2 / pow( ShearM, 2) * Effect * EffectDeriv * std::sin(3.0*LodeAngle);

         //C3 = -3 * sqrt(3.0) * Effect * EffectDeriv / ( pow( ShearM, 2.0) * std::cos( 3.0*LodeAngle) ); // removing cos
         C3 = -3 * sqrt(3.0) * Effect * EffectDeriv / ( pow( ShearM, 2) ); 
         C3 /= J2;
      }
      Vector ThisDerivative =  C2* V2 + C3*V3;

      /*std::cout << " LODE " << LodeAngle <<" LODE " << LodeAngle * 180.0 / GetPI() << " EFFECT " << Effect <<  " DERIVATIVE " << EffectDeriv * std::cos( 3.0 * LodeAngle)  << std::endl;
        std::cout << " PREVIOUS ANAL DERIVATIVE " << rYieldFunctionD << std::endl;
        std::cout << " AND THIS NEW C2 " << C2 << " and C3 " << C3 << std::endl;
        std::cout << " V2 " << V2 << " V3 " << V3 << std::endl;*/
      rYieldFunctionD += ThisDerivative;


   }

   //************************* smoothingInvariants of something ***********
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

      rB = -1.0 * ( Sign* std::sin(SmoothingAngle) + std::sin(FrictionAngle)*std::cos(SmoothingAngle) / sqrt(3.0) ) / ( 3.0*std::cos(3.0*SmoothingAngle) );



   }
   double CamClayYieldCriterion::GetSmoothingLodeAngle()
   {
      return 27.0*GetPI()/180.0;
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
