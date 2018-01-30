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
#include "includes/ublas_interface.h"

#include "custom_constitutive/custom_flow_rules/cam_clay_explicit_plastic_flow_rule.hpp"

#include "pfem_solid_mechanics_application_variables.h"

// TO BE REMOVED. This is a special case of Borja


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

      /*    std::cout << " COMPUTING STRESS " << rHenckyStrainVector << std::endl;
            std::cout << " DeviatoricStrainVector " << DeviatoricStrainVector << std::endl;
            std::cout << " Dev Stress " << rKirchhoffStressVector << std::endl;
            std::cout << " VOL STRAIN " << VolumetricStrain << std::endl;
            std::cout << " VOL STRESS " << MeanStress << std::endl;
            std::cout << " " << std::endl;*/
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

      double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
      double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];

      double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
      double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
      ReferencePressure /= OCR;    

      double DeviatoricStrain2Norm = 0.0;
      for (unsigned int i = 0; i < 3; ++i)
         DeviatoricStrain2Norm += pow(rDeviatoricStrainVector(i), 2);

      for (unsigned int i = 3; i < 6; ++i)
         DeviatoricStrain2Norm += 2.0*pow(rDeviatoricStrainVector(i)/2.0, 2);

      rMeanStress = -ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope) * (1.0 + 1.0*AlphaShear*DeviatoricStrain2Norm/SwellingSlope);


      /*   std::cout << " -------" << std::endl;
           std::cout << " VOLUMETRICMOD " << -ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope) * ( 1.0 + AlphaShear*DeviatoricStrain2Norm / SwellingSlope) << std::endl;
           std::cout << "      ShearNor " << DeviatoricStrain2Norm << std::endl;
           std::cout << "      ShearTer " << AlphaShear*DeviatoricStrain2Norm/SwellingSlope << std::endl;
           std::cout << "      VolStrain" << rVolumetricStrain << std::endl;
           std::cout << "      MeanStres" << rMeanStress << std::endl; 
           std::cout << " -------" << std::endl;
       */

   }




   void CamClayExplicitFlowRule::CalculateDeviatoricStress(const double& rVolumetricStrain, const Vector & rDeviatoricStrainVector, Vector& rDeviatoricStress)
   {
      double ReferencePressure = mpYieldCriterion->GetHardeningLaw().GetProperties()[PRE_CONSOLIDATION_STRESS];
      double OCR = mpYieldCriterion->GetHardeningLaw().GetProperties()[OVER_CONSOLIDATION_RATIO];
      ReferencePressure /= OCR;    
      double SwellingSlope = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
      double AlphaShear = mpYieldCriterion->GetHardeningLaw().GetProperties()[ALPHA_SHEAR];


      rDeviatoricStress = rDeviatoricStrainVector;
      rDeviatoricStress *= 2.0*AlphaShear*ReferencePressure* std::exp(-rVolumetricStrain / SwellingSlope);

      for (unsigned int i = 3; i<6; ++i)
         rDeviatoricStress(i) /= 2.0;  // BECAUSE VOIGT NOTATION

      /*std::cout << " SHEAR MODULUS " << 2.0*AlphaShear*ReferencePressure*std::exp( -rVolumetricStrain / SwellingSlope) << std::endl;
        std::cout << "     DevStrain " << rDeviatoricStrainVector << std::endl;
        std::cout << "     VolumeStra" << rVolumetricStrain << std::endl;
        std::cout << "     rDeviatStr" << rDeviatoricStress << std::endl;
       */

   }


   void CamClayExplicitFlowRule::ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix )
   {


      Matrix FourthOrderIdentity = ZeroMatrix(6,6);
      for (unsigned int i = 0; i<3; ++i)
         FourthOrderIdentity(i,i) = 1.0;

      for (unsigned int i = 3; i<6; ++i)
         FourthOrderIdentity(i,i) = 0.50;
      // VOIGT NOTATION AND NOT KELVIN

      Matrix IdentityCross = ZeroMatrix(6,6);
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


      rElasticMatrix  = (-1.0/SwellingSlope)*MeanStress*IdentityCross;
      rElasticMatrix += 2.0*AlphaShear*ReferencePressure*std::exp(-VolumetricStrain/SwellingSlope)*(FourthOrderIdentity - (1.0/3.0)*IdentityCross);


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

      //   if ( fabs(VolumetricStrain) > 0.50) {
      //     for (unsigned int i = 0; i<6; ++i) {
      //       rElasticMatrix(i,i) += 0000.1;
      //     }
      //   }
      //std::cout << "ELASTIC MATRIX " << rElasticMatrix << "ElasticVector " << rElasticStrainVector << " StressV " << StressVector << "Mean Stress " << MeanStress << std::endl;

   }


   void CamClayExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
   {

      double MeanStress;
      this->CalculateMeanStress(rHenckyStrainVector, MeanStress);

      double PreconsolidationStress = mpYieldCriterion->GetHardeningLaw().CalculateHardening(PreconsolidationStress, rAlpha);
      double SwellingSlope          = mpYieldCriterion->GetHardeningLaw().GetProperties()[SWELLING_SLOPE];
      double OtherSlope             = mpYieldCriterion->GetHardeningLaw().GetProperties()[NORMAL_COMPRESSION_SLOPE];


      rH = (2.0*MeanStress-PreconsolidationStress) ;
      rH *= (-MeanStress);
      rH *= PreconsolidationStress/ ( OtherSlope - SwellingSlope);

   }


   void CamClayExplicitFlowRule::CalculatePlasticPotentialDerivatives( const Vector& rStressVector, Vector& rFirstDerivative, Matrix & rSecondDerivative)
   {

      rFirstDerivative = ZeroVector(1);
      rSecondDerivative = ZeroMatrix(1,1);
      double M = mpYieldCriterion->GetHardeningLaw().GetProperties()[CRITICAL_STATE_LINE] ;
      double Friction = mpYieldCriterion->GetHardeningLaw().GetProperties()[INTERNAL_FRICTION_ANGLE];


      if ( fabs( Friction) > 1e-5)
         return;

      Vector StressV = rStressVector;
      double MeanP = 0;
      for (unsigned int i=0 ; i < 3; i++)
         MeanP += StressV(i) / 3.0;

      for (unsigned int i = 0; i < 3; i++)
         StressV(i) -= MeanP;

      // COMPUTE J2;
      double J2InvSQ = 0;
      for (unsigned int i = 0; i < 3; i++)
         J2InvSQ += pow( StressV(i), 2);

      for (unsigned int i = 0; i < 3; i++)
         J2InvSQ += 2.0 * pow( StressV(i), 2);

      J2InvSQ = sqrt( J2InvSQ/2.0);

      if ( J2InvSQ < 1E-5)
         return;

      // J2 firstDeriv
      Vector C2Vector = StressV;
      for (unsigned int i = 3; i < 6; i++)
         C2Vector(i) *= 2.0;
      C2Vector /= 2.0* J2InvSQ; 

      // C2Matrix
      Matrix C2Matrix = ZeroMatrix(6,6);
      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            double times = 0.5;
            if ( i > 2)
               times *= 2.0;
            if ( j > 2 )
               times *= 2.0;
            C2Matrix(i,j) -= times* StressV(i)*StressV(j) / pow( J2InvSQ, 2);
         }
      }

      for (unsigned int i = 0; i < 3; i++) {
         C2Matrix(i,i) += 1.0;
         for (unsigned int j = 0; j < 3; j++)  {
            C2Matrix(i,j) -= 1.0/3.0;
         }
      }
      for (unsigned int i = 3; i < 6; i++)
         C2Matrix(i,i) += 2.0;

      C2Matrix /= 2.0 * J2InvSQ;

      rSecondDerivative = 6.0 / ( M*M) * J2InvSQ * C2Matrix;

      Vector C1Der = ZeroVector(6);
      for (unsigned int i = 0; i < 3; i++)
         C1Der(i) = 1.0/3.0;

      for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
            rSecondDerivative(i,j) += 6.0/ (M*M) * C2Vector(i)*C2Vector(j) + 2.0*C1Der(i)*C1Der(j);
         }
      }

   }

   void CamClayExplicitFlowRule::save( Serializer& rSerializer) const 
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
   }

   void CamClayExplicitFlowRule::load( Serializer& rSerializer)
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )

   }


} //end namespace kratos
