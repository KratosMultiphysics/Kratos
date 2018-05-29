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
#include "custom_constitutive/custom_yield_criteria/new_tresca_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


   //*******************************CONSTRUCTOR******************************************
   //************************************************************************************
   NewTrescaYieldCriterion::NewTrescaYieldCriterion()
      :TrescaYieldCriterion()
   {

   }

   //*****************************INITIALIZATION CONSTRUCTOR*****************************
   //************************************************************************************

   NewTrescaYieldCriterion::NewTrescaYieldCriterion(HardeningLawPointer pHardeningLaw)
      :TrescaYieldCriterion(pHardeningLaw)
   {

   }

   //*******************************ASSIGMENT OPERATOR***********************************
   //************************************************************************************

   NewTrescaYieldCriterion& NewTrescaYieldCriterion::operator=(NewTrescaYieldCriterion const& rOther)
   {
      TrescaYieldCriterion::operator=(rOther);
      return *this;
   }

   //*******************************COPY CONSTRUCTOR*************************************
   //************************************************************************************

   NewTrescaYieldCriterion::NewTrescaYieldCriterion(NewTrescaYieldCriterion const& rOther)
      :TrescaYieldCriterion(rOther)
   {

   }


   //********************************DESTRUCTOR******************************************
   //************************************************************************************

   NewTrescaYieldCriterion::~NewTrescaYieldCriterion()
   {
   }



   //************************* CALCULATE YIELD FUNCTION  ******************
   //**********************************************************************

   double& NewTrescaYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha)
   {

      KRATOS_TRY

      double Su = this->GetHardeningLaw().GetProperties()[YIELD_STRESS];
      double AlphaS = this->GetHardeningLaw().GetProperties()[LAMBDA]; // quotient between the compression and extension undrained shear strenght
      if ( AlphaS < 1.0 || AlphaS > 2.0)
         AlphaS = 1.0;
      double smoothing = 10.0; // to smooth the sharp corners


      Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( rStressVector); 

      TrescaStressInvariants StressInvariants;
      this->CalculateStressInvariants( rStressVector, StressInvariants);

      double J2 = StressInvariants.J2InvSQ;
      J2 = pow(J2, 2);
      double lode = -StressInvariants.LodeAngle;

      double   t2 = sqrt(3.0); 
      double   t3 = sin(lode); 
      double   t4 = AlphaS*3.0; 
      double   t5 = t4-3.0; 
      double   t6 = t3*t5; 
      double   t7 = cos(lode); 
      double   t8 = AlphaS+1.0; 
      double   t9 = t2*t7*t8; 
      double   t10 = 3.141592653589793*(2.0/3.0); 
      double   t11 = lode+t10; 
      double   t13 = 3.141592653589793*(1.0/3.0); 
      double   t12 = lode-t13; 
      double   t14 = lode+t13; 
      double   t15 = 3.141592653589793*(4.0/3.0); 
      double   t16 = lode+t15; 
      double   t0 = -Su+(sqrt(J2)*t2*log(exp(smoothing*(t5*sin(t11)+t2*t8*cos(t11)))+exp(-smoothing*(t5*sin(t12)-t2*t8*cos(t12)))+exp(-smoothing*(t5*sin(t14)-t2*t8*cos(t14)))+exp(smoothing*(t5*sin(t16)+t2*t8*cos(t16)))+exp(smoothing*(t6+t9))+exp(smoothing*(t6-t9)))*(1.0/6.0))/smoothing; 
      rStateFunction = t0; 


      return rStateFunction; 

      KRATOS_CATCH("")
   }




   //************************* YIELD FUNCTION DERIVATIVE ******************
   //**********************************************************************
   void NewTrescaYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
   {
      KRATOS_TRY

      double pertur = pow(2, -25);

      rYieldFunctionD = ZeroVector(6);

      Vector Stress1 = ZeroVector(6);
      Vector Stress2 = ZeroVector(6);
      
      double F1, F2;
      
      for (unsigned int i = 0; i < 6; i++) {
         Stress1 = rStressVector;
         Stress1(i) += pertur;
         F1 = this->CalculateYieldCondition( F1, Stress1, rAlpha);
         Stress2 = rStressVector;
         Stress2(i) -= pertur;
         F2 = this->CalculateYieldCondition( F2, Stress2, rAlpha);
         double denominador = Stress1(i)-Stress2(i);

         rYieldFunctionD(i) = (F1-F2)/(denominador);
      }

      for (unsigned int i = 3; i <6; i++)
         rYieldFunctionD(i) *= 2.0;

      KRATOS_CATCH("")

   }
   void NewTrescaYieldCriterion::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, TrescaYieldCriterion )
   }

   void NewTrescaYieldCriterion::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, TrescaYieldCriterion )
   }


}  // namespace Kratos.
