//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>
#include<cmath>

// External includes
#include "includes/ublas_interface.h"

// Project includes
#include "custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "pfem_solid_mechanics_application_variables.h"


#include "custom_utilities/stress_invariants_utilities.hpp"
namespace Kratos
{


   // ******************* DEFAULT CONSTRUCTOR ***************************
   NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule()
      :FlowRule()
   {

   }

   //*****************************INITIALIZATION CONSTRUCTOR*****************************
   //************************************************************************************
   NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
      :FlowRule(pYieldCriterion)
   {

   }

   // ******************* COPY CONSTRUCTOR  *******************************************
   NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule(const NonAssociativeExplicitPlasticFlowRule& rOther)
      :FlowRule(rOther),
       mPlasticMultiplierVelocity(rOther.mPlasticMultiplierVelocity)
   {

   }

   NonAssociativeExplicitPlasticFlowRule& NonAssociativeExplicitPlasticFlowRule::operator=( NonAssociativeExplicitPlasticFlowRule const & rOther)
   {
      FlowRule::operator=(rOther);
      mPlasticMultiplierVelocity = rOther.mPlasticMultiplierVelocity; 
      return *this;
   }

   FlowRule::Pointer NonAssociativeExplicitPlasticFlowRule::Clone() const
   {
      FlowRule::Pointer p_clone(new NonAssociativeExplicitPlasticFlowRule(*this));
      return p_clone;
   }


   NonAssociativeExplicitPlasticFlowRule::~NonAssociativeExplicitPlasticFlowRule()
   {
   }


   // ******* Initilialize Material (due to the new variable) ******
   void NonAssociativeExplicitPlasticFlowRule::InitializeMaterial( YieldCriterionPointer & pYieldCriterion, HardeningLawPointer & pHardeningLaw, const Properties & rMaterialProperties)
   {
      mPlasticMultiplierVelocity = 0.0;
      FlowRule::InitializeMaterial( pYieldCriterion, pHardeningLaw, rMaterialProperties);
   }

   void NonAssociativeExplicitPlasticFlowRule::InitializeMaterial( const Properties & rMaterialProperties)
   {
      mPlasticMultiplierVelocity = 0.0;
      FlowRule::InitializeMaterial( rMaterialProperties);
   }
   



   bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
   {
      rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;

      bool PlasticityActive = false;

      int StressIntTechnique = 0; // mpYieldCriterion->GetHardeningLaw().GetProperties()[STRESS_INTEGRATION_STRATEGY];

      if ( ! rReturnMappingVariables.Options.Is( IMPLEX_ACTIVE) )
      {
         if ( StressIntTechnique == 2)
         {
            // RETURN MAPPING ALGORITHM
            PlasticityActive = CalculateReturnMappingImpl( rReturnMappingVariables, rIncrementalDeformationGradient, rStressMatrix, rNewElasticLeftCauchyGreen); 
         }
         else {
            // EXPLICIT STRESS INTEGRATION WITH SUBSTEPPING
            PlasticityActive = CalculateReturnMappingExpl( rReturnMappingVariables, rIncrementalDeformationGradient, rStressMatrix, rNewElasticLeftCauchyGreen);
         }
         //return PlasticityActive;
      }
      if ( rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE) ) {
         if ( StressIntTechnique == 4)
         {
            // IMPLEX STRESS INTEGRATION
            PlasticityActive = CalculateReturnMappingImplex2( rReturnMappingVariables, rIncrementalDeformationGradient, rStressMatrix, rNewElasticLeftCauchyGreen);
         }
         else {
            // MODIFIED (explicit) IMPLEX STRESS INTEGRATION
            PlasticityActive = CalculateReturnMappingImplex( rReturnMappingVariables, rIncrementalDeformationGradient, rStressMatrix, rNewElasticLeftCauchyGreen);
         }
         return PlasticityActive;
      }


      if ( PlasticityActive == true && false )
      {
         std::cout << " STUPID FUNCTION TO EVALUATE THE YIELD FUNCTION DERIVATIVE NUMERICALLY AND ANALOGICALLY " << std::endl;
         Vector StressVector;
         StressVector = MathUtils<double>::StressTensorToVector( rStressMatrix, 6);

         // GET ALPHA
         InternalVariables PlasticVariables = mInternalVariables; 
         double Alpha = PlasticVariables.EquivalentPlasticStrain; 

         // ANALYTHICAL DERIVATIVE
         Vector YieldFunctionD; 
         mpYieldCriterion->CalculateYieldFunctionDerivative( StressVector, YieldFunctionD, Alpha);

         std::cout << " ANALYTHICAL DERIVATIVE " << std::endl;
         std::cout << YieldFunctionD << std::endl;

         // NUMERICAL DERIVATIVE
         double F0; 
         F0 = mpYieldCriterion->CalculateYieldCondition( F0, StressVector, Alpha);
         double delta = 5.0e-7;
         double F1;
         for (unsigned int i = 0; i < 6; i++)
         {
            StressVector(i) += delta; 
            F1 = mpYieldCriterion->CalculateYieldCondition(F1, StressVector, Alpha);
            YieldFunctionD(i) = (F1-F0) / delta; 
            StressVector(i) -= delta; 
         }
         std::cout << " NUMERICAL DERIVATYIVE"  << std::endl;   
         std::cout << YieldFunctionD << std::endl;
         std::cout << std::endl;

         // I MOVE THE STRESS STATE
         /*double MeanP = 0; 
           for ( int i = 0; i < 3; i++)
           MeanP += StressVector(i);
           for (int i = 0; i < 3; i++)
           StressVector(i) += ( 40 - MeanP / 3 );*/

         Matrix EigenVectors;
         Vector EigenValuesS; 
         SolidMechanicsMathUtilities<double>::EigenVectors( rStressMatrix, EigenVectors, EigenValuesS);

         StressVector = ZeroVector(6);
         for (unsigned int i = 0; i < 3; i++)
            StressVector(i) = EigenValuesS(i);
         mpYieldCriterion->CalculateYieldFunctionDerivative( StressVector, YieldFunctionD, Alpha);

         /*std::cout << " EXAMPLE " << EigenValuesS(0) << " " << EigenValuesS(1) << " " << EigenValuesS(2) ;

           for (unsigned int i = 0; i < 6; i++)
           std::cout << " " << YieldFunctionD(i) ; 
           std::cout << std::endl;*/

         Vector ThisStress = ZeroVector(6);
         ThisStress(0) = -60;
         ThisStress(2) = -20;

         for (unsigned int i = 0; i < 1000; i++)
         {
            ThisStress(1) = ThisStress(0) - (ThisStress(0) - ThisStress(2) )* double(i) / (1000.0 - 1.0);
            mpYieldCriterion->CalculateYieldFunctionDerivative( ThisStress, YieldFunctionD, Alpha);
            std::cout << " EXAMPLE " << ThisStress(0) << " " << ThisStress(1) << " " << ThisStress(2) ;

            for (unsigned int i = 0; i < 6; i++)
               std::cout << " " << YieldFunctionD(i) ; 
            std::cout << std::endl;
            double I1, I2, I3;
            StressInvariantsUtilities::CalculateStressInvariants( ThisStress, I1, I2, I3);
            std::cout << " THIS INVARIANTS " << I1 << " " << I2 << " " << I3 << std::endl;

            // try to get the numerical derivative
            double F0; 
            F0 = mpYieldCriterion->CalculateYieldCondition( F0, ThisStress, Alpha);
            double delta = 5.0e-7;
            double F1;
            for (unsigned int i = 0; i < 6; i++)
            {
               ThisStress(i) += delta; 
               F1 = mpYieldCriterion->CalculateYieldCondition(F1, ThisStress, Alpha);
               YieldFunctionD(i) = (F1-F0) / delta; 
               ThisStress(i) -= delta; 
            }
            std::cout << " NUMERICAL DERIVATYIVE"  << std::endl;   
            std::cout << YieldFunctionD << std::endl;
            std::cout << std::endl;

         }



         KRATOS_THROW_ERROR(std::logic_error, " I HATE YOUUU, its a joke " , "")


      }

      return PlasticityActive;

   }


   // MODIFIED (explicit) IMPLEX STRESS INTEGRATION
   bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMappingImplex(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
   {

      bool PlasticityActive = false;
      InternalVariables PlasticVariables = mInternalVariables;
      rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;
      rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain; 
      double PlasticMultiplier = mPlasticMultiplierVelocity;

      Matrix PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
      Vector PreviousStressVector;
      this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, PreviousStressVector);

      Vector NewStressVector;

      double Tolerance = 1E-6;

      // STEP 1. PRODUCE THE EXPLICIT

      // Actualize the Plastic Multiplier (obs: we save Plastic*Delta time, this way, if it changes, ...)
      PlasticMultiplier *= rReturnMappingVariables.DeltaTime;

      // TRIAL ONLY ELASTIC STAGE
      Matrix UpdateMatrix = rDeltaDeformationGradient; 
      UpdateMatrix = rDeltaDeformationGradient;
      rNewElasticLeftCauchyGreen = prod(UpdateMatrix, PreviousElasticLeftCauchyGreen);
      rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix));

      this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, NewStressVector);
      rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, NewStressVector, rReturnMappingVariables.DeltaGamma);

      // IF NOT ELASTIC, THERE IS CLEARLY PLASTIC DEFORMATION .with explicit value.
      if ( rReturnMappingVariables.TrialStateFunction > Tolerance)
      {
         if ( fabs(PlasticMultiplier) < 1e-12)
            PlasticMultiplier = 1e-11;

         AuxiliarDerivativesStructure  AuxiliarDerivatives;
         Vector HenckyStrain = ConvertCauchyGreenTensorToHenckyStrain (PreviousElasticLeftCauchyGreen );

         this->UpdateDerivatives( HenckyStrain, AuxiliarDerivatives, rReturnMappingVariables.DeltaGamma );
         Vector AuxiliarVector =  -PlasticMultiplier * AuxiliarDerivatives.PlasticPotentialD / 2.0;

         UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor( AuxiliarVector);

         UpdateMatrix = prod( rDeltaDeformationGradient, UpdateMatrix);
         rNewElasticLeftCauchyGreen = prod( UpdateMatrix, PreviousElasticLeftCauchyGreen);
         rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix) );

         // ACTUALIZE HARDENING PARAMENTERS (Volumetric) 
         double AlphaUpdate = 0;
         for (unsigned int i = 0; i < 3; i++) {
            rReturnMappingVariables.DeltaGamma += PlasticMultiplier * AuxiliarDerivatives.PlasticPotentialD(i);
            AlphaUpdate += PlasticMultiplier * AuxiliarDerivatives.PlasticPotentialD(i);
         }

         double aux = 0.0;
         for (unsigned int j = 0; j < 3; ++j)
            aux += pow(PlasticMultiplier*AuxiliarDerivatives.PlasticPotentialD(j) - AlphaUpdate/3.0, 2);

         for (unsigned int j = 3; j < 6; ++j)
            aux  += 2.0*pow(PlasticMultiplier*AuxiliarDerivatives.PlasticPotentialD(j) / 2.0, 2);
         rReturnMappingVariables.IncrementalPlasticShearStrain += sqrt( aux);

         PlasticityActive = true;
      }

      if ( PlasticityActive )
      {
         rReturnMappingVariables.NormIsochoricStress = 0.0*PlasticMultiplier; // no, perque hi guardo el multiplicador plastic per poder fer la matriu tangent consistent.
      }
      else {
         rReturnMappingVariables.NormIsochoricStress = 0.0;
      }
      rStressMatrix = ComputeKirchhoffStressMatrix( rNewElasticLeftCauchyGreen);
      rReturnMappingVariables.Options.Set(PLASTIC_REGION, PlasticityActive);

      /*std::cout << " PREVIOUS STRESS " << PreviousStressVector << std::endl;
        std::cout << " STATE FUNCTION " << rReturnMappingVariables.TrialStateFunction << std::endl;
        std::cout << " PREIVIOUS UPDATE MATRIX " << rDeltaDeformationGradient << std::endl; 
        std::cout << " PLASTIC ACTIVE " << PlasticityActive << std::endl;
        std::cout << " UPDATE MATRIX " << UpdateMatrix << std::endl;
        std::cout << " PLASTIC MULTIPLIER " << PlasticMultiplier << " and R " << rReturnMappingVariables.DeltaTime << std::endl;
        std::cout << " PREVIOUS STRAIN " << PreviousElasticLeftCauchyGreen << std::endl;
        std::cout << " STRAIN " << rNewElasticLeftCauchyGreen << std::endl;
        std::cout << " HARD LIKE " << rReturnMappingVariables.DeltaGamma << std::endl;
        std::cout << " EXIT STRESS " << rStressMatrix << std::endl;
        std::cout << " NEW THAT " << rReturnMappingVariables.NormIsochoricStress << std::endl;
        std::cout << " " << std::endl;*/
      return PlasticityActive;


   }

   // IMPLEX ALGORITHM
   bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMappingImplex2(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
   {

      bool PlasticityActive = false;
      InternalVariables PlasticVariables = mInternalVariables;
      rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;
      rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain; 
      double PlasticMultiplier = mPlasticMultiplierVelocity;

      Matrix PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
      Vector PreviousStressVector;

      Vector NewStressVector;

      //double Tolerance = 1E-6;

      // 2- Check for yield Condition (at the beggining)
      this->CalculateKirchhoffStressVector( PreviousElasticLeftCauchyGreen, PreviousStressVector);
      rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, rReturnMappingVariables.DeltaGamma );


      // STEP 1. PRODUCE THE EXPLICIT
      PlasticMultiplier *= rReturnMappingVariables.DeltaTime;

      // TRIAL ONLY ELASTIC STAGE
      Matrix UpdateMatrix = rDeltaDeformationGradient; 
      rNewElasticLeftCauchyGreen = prod(UpdateMatrix, PreviousElasticLeftCauchyGreen);
      rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix));

      this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, NewStressVector);
      rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, NewStressVector, rReturnMappingVariables.DeltaGamma);

      // THERE WAS A PLASTIC STAGE
      if (  ( (PlasticMultiplier) > 1e-11) ) //&& rReturnMappingVariables.TrialStateFunction > 0 )
      {

         // Define some variables previous to iterate
         AuxiliarDerivativesStructure  AuxiliarDerivatives;
         Vector HenckyStrainTrial = ConvertCauchyGreenTensorToHenckyStrain ( rNewElasticLeftCauchyGreen );
         Vector DeltaHencky; 
         Vector HenckyStrain = HenckyStrainTrial; // this is the unkwonw
         Vector Residual; 
         Matrix SystemMatrix, InverseSystemMatrix, ElasticMatrix; 
         double ErrorNorm;
         bool is_converged = false; 


         // EXPLICIT PREDICTION OF THE "hardening" PARAMTER (volumetric)
         Vector HenckyStrainPrev = ConvertCauchyGreenTensorToHenckyStrain( PreviousElasticLeftCauchyGreen);
         this->UpdateDerivatives( HenckyStrainPrev, AuxiliarDerivatives, rReturnMappingVariables.DeltaGamma );
         for (unsigned int i = 0; i < 3; i++) {
            rReturnMappingVariables.DeltaGamma += PlasticMultiplier * AuxiliarDerivatives.PlasticPotentialD(i);
         }


         // BEGINING OF SOLVING THE NON-LINEAR THING 
         this->UpdateDerivatives( HenckyStrain, AuxiliarDerivatives, rReturnMappingVariables.DeltaGamma );

         for (unsigned int iter = 0; iter < 150; iter ++)
         {

            this->UpdateDerivatives( HenckyStrain, AuxiliarDerivatives, rReturnMappingVariables.DeltaGamma );
            Residual = HenckyStrain + PlasticMultiplier * AuxiliarDerivatives.PlasticPotentialD - HenckyStrainTrial;

            ErrorNorm = 0; 
            for ( int k = 0; k < 6; ++k)
            {
               ErrorNorm += pow( Residual(k), 2);
            }
            ErrorNorm = sqrt( ErrorNorm);


            if (  ErrorNorm == ErrorNorm +1.0 )
            {
               is_converged = false; 
               std::cout << " NAN CORRECTION2 " << std::endl;
               break;
            }
            else if ( std::isnan( ErrorNorm) )
            {
               is_converged = false; 
               std::cout << " NAN CORRECTION " << std::endl;
               break;
            } 

            if ( ( ErrorNorm < 5e-5) && iter > 0)
            {
               is_converged = true;
               break; 
            }


            SystemMatrix = ZeroMatrix(6,6);
            for (int ii = 0; ii < 6; ii++)
            {
               SystemMatrix(ii,ii) = 1.0;
            }

            this->ComputeElasticMatrix(HenckyStrain, ElasticMatrix);
            Matrix ThisUpdate = PlasticMultiplier * prod( AuxiliarDerivatives.PlasticPotentialDD, ElasticMatrix);
            SystemMatrix += ThisUpdate; 

            InverseSystemMatrix = ZeroMatrix(6,6);
            SolidMechanicsMathUtilities<double>::InvertMatrix(SystemMatrix, InverseSystemMatrix);

            DeltaHencky = prod( InverseSystemMatrix, Residual);
            Vector HenckyStrainUpdate = -DeltaHencky;

            HenckyStrain += HenckyStrainUpdate;

         }

         if ( is_converged == false )
         {
            std::cout << " trial " ;
            this->UpdateDerivatives( HenckyStrainTrial, AuxiliarDerivatives, rReturnMappingVariables.DeltaGamma );

            HenckyStrain = HenckyStrainTrial - PlasticMultiplier * AuxiliarDerivatives.PlasticPotentialD;
         }

         double aux = 0.0, AlphaUpdate = 0.0;
         for (unsigned int j = 0; j < 3; ++j)
            aux += pow(PlasticMultiplier*AuxiliarDerivatives.PlasticPotentialD(j) - AlphaUpdate/3.0, 2);

         for (unsigned int j = 3; j < 6; ++j)
            aux  += 2.0*pow(PlasticMultiplier*AuxiliarDerivatives.PlasticPotentialD(j) / 2.0, 2);
         rReturnMappingVariables.IncrementalPlasticShearStrain += sqrt( aux);

         rNewElasticLeftCauchyGreen = ConvertHenckyStrainToCauchyGreenTensor( HenckyStrain ); 

         PlasticityActive = true;

      }

      if ( PlasticityActive )
      {
         rReturnMappingVariables.NormIsochoricStress = PlasticMultiplier;
      }
      else {
         rReturnMappingVariables.NormIsochoricStress = 0.0;
      }
      rStressMatrix = ComputeKirchhoffStressMatrix( rNewElasticLeftCauchyGreen);
      rReturnMappingVariables.Options.Set(PLASTIC_REGION, PlasticityActive);


      return PlasticityActive;


   }

   // RETURNS A STRESS POINT TO THE YIELD SURFACE ACCORDING TO THE RETURN MAPPING ALGORITHMS (Some terms of the system matrix for the MCC are missing, though, thay may be not very revelant)
   void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface4( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rDrift, const double& rTolerance)
   {

      AuxiliarDerivativesStructure  AuxiliarDerivatives;
      Matrix ElasticMatrix;
      double H;
      Matrix  UpdateMatrix;
      Vector StressVector;
      double AlphaTrial = rReturnMappingVariables.DeltaGamma; 
      double Alpha = AlphaTrial; 
      rReturnMappingVariables.DeltaBeta = 0.0;  // THIS STEP PLASTIC MULTIPLIER

      Vector HenckyTrialElastic = ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);
      Vector HenckyElastic = HenckyTrialElastic; 


      Vector R1 ;
      double R3 = 0, R2 = 0;
      double Gamma = 0.0, ResidualNorm = 0, PreviousResidualNorm = 1000000000; 
      Vector Residual = ZeroVector(8);
      Matrix SystemMatrix, InverseSystemMatrix, AuxMatrix; 
      Vector AuxVector, deltaX = ZeroVector(8);  


      if ( false )
      {
         this->UpdateDerivatives( HenckyElastic, AuxiliarDerivatives, Alpha);

         std::cout << " CHECK THE DERIVATIVE " << std::endl;
         Vector KirS ;
         this->CalculateKirchhoffStressVector( HenckyElastic, KirS);
         Vector StressV = KirS;
         Vector NumDer = ZeroVector(6);
         double delta = 1e-7;
         double Ryield; 
         Ryield = mpYieldCriterion->CalculateYieldCondition( Ryield, StressV, Alpha);

         for (unsigned int m = 0; m < 6; m++) {
            double yield ;
            StressV = KirS;
            StressV(m) += delta; 
            yield = mpYieldCriterion->CalculateYieldCondition( yield, StressV, Alpha);
            NumDer(m) = (yield-Ryield) / delta;
         }
         std::cout << " ANALYTICAL " << AuxiliarDerivatives.YieldFunctionD << std::endl;
         std::cout << " NUMERICAL  " << NumDer << std::endl;
         std::cout << " " << std::endl;
         std::cout << " " << std::endl;
         // AND NOW THE SECOND DERIVATIVE.....
         delta = 1e-5;
         for (unsigned int m = 0; m < 6; m++)
         {
            Vector Der;
            StressV = KirS;
            StressV(m) += delta; 
            mpYieldCriterion->CalculateYieldFunctionDerivative( StressV, Der, Alpha);
            for (unsigned int kk = 0; kk < 6; kk++)
               Der(kk) = (Der(kk) - AuxiliarDerivatives.YieldFunctionD(kk) ) / delta;

            std::cout << " comp " << m ;
            for ( int tt = 0; tt < 6; tt ++)
               std::cout << " " << Der(tt) ;
            std::cout << std::endl;

         }
         std::cout << " ANALYTICAL " << std::endl;;
         for (unsigned t = 0; t < 6; t++) {
            for (unsigned q = 0; q < 6; q++) {
               std::cout << " "  <<AuxiliarDerivatives.PlasticPotentialDD(t, q);
            }
            std::cout << " " << std::endl;
         }
         std::cout << " ENDSECONDDERIVATIVE " << std::endl;
         std::cout << " " << std::endl;
         std::cout << " " << std::endl;
         std::cout << " " << std::endl;
         std::cout << " " << std::endl;


      }


      // ITERATIONS TO PERFORM THE RETURN MAPPING

      for (unsigned int iter = 0; iter < 150; ++iter) {

         // ADD SOME SORT OF
         for (unsigned int i = 0; i < 6; i++)
            HenckyElastic(i) -= deltaX(i);

         Alpha -= deltaX(6);
         Gamma -= deltaX(7);

         this->UpdateDerivatives( HenckyElastic, AuxiliarDerivatives, Alpha);

         // COMPUTE THE RESIDUAL;
         R1 = HenckyElastic + Gamma * AuxiliarDerivatives.PlasticPotentialD - HenckyTrialElastic; 

         double AuxP = 0;
         for (unsigned int i = 0; i < 3; i++)
            AuxP += AuxiliarDerivatives.PlasticPotentialD(i);
         R2 = Alpha - AlphaTrial - Gamma *  AuxP; 

         this->CalculateKirchhoffStressVector( HenckyElastic , rStressVector);
         R3 = mpYieldCriterion->CalculateYieldCondition( R3, rStressVector, Alpha);

         // PUT IT TOGUETHER and COMPUTE NORM
         for (unsigned int i = 0; i < 6; i++)
            Residual(i) = R1(i);
         Residual(6) = R2;
         Residual(7) = R3; 

         ResidualNorm = 0;
         for (unsigned int i = 0; i < 8; i++)
            ResidualNorm += pow( Residual(i), 2);
         ResidualNorm = sqrt(ResidualNorm); 

         //std::cout << "   RRMM ITER " << iter << " NORM " << ResidualNorm << std::endl;;
         if ( ( ResidualNorm > PreviousResidualNorm) || true )
         {
            //std::cout << " YOU HAVE TO PERFORM SOME SORT OF SOMETHING " << std::endl;
            //std::cout << " " << ResidualNorm << " " << PreviousResidualNorm << std::endl;
            //std::cout << " iteration " << iter << std::endl;
            PerformSomeSortOfLineSearch( HenckyElastic, Gamma, Alpha, HenckyTrialElastic, AlphaTrial, ResidualNorm, PreviousResidualNorm, deltaX);
            //std::cout << " " << ResidualNorm << " " << PreviousResidualNorm << std::endl;
            this->UpdateDerivatives( HenckyElastic, AuxiliarDerivatives, Alpha);
         }
         //std::cout << "   RRMM ITER " << iter << " NORM " << ResidualNorm << std::endl;;

         if ( iter > 140)
         { // WRITTING THINGS IN THAT

            std::cout << " " << iter << "  " << ResidualNorm ;

            this->CalculateKirchhoffStressVector( HenckyElastic , rStressVector);
            Matrix SM = MathUtils<double>::StressVectorToTensor( rStressVector);
            double mP = 0;
            for (int i = 0; i < 3; i++)
               mP += SM(i,i) / 3.0;
            for (int i = 0; i < 3; i++)
               SM(i,i) -= mP;
            double J2 = 0;
            for (int i = 0; i < 3; i++) {
               for (int j = 0; j < 3; j++) {
                  J2 += SM(i,j)*SM(i,j) / 2;
               }
            }
            J2 = sqrt(J2);

            double Lode = MathUtils<double>::Det(SM);
            Lode *= -3.0*sqrt(3.0) / 2.0 / pow( J2, 3);
            Lode = std::asin(Lode) / 3.0;
            Lode *= 180.0 / 3.141593;
            std::cout << " " << mP << " " << J2 << " " << Lode << std::endl;


         }



         if ( ResidualNorm < 1E-5)
         {
            break; 
         }


         // COMPUTE THE MATRIX OF THE SYSTEM ( falten termes que els tinc en una llibreta)
         SystemMatrix = ZeroMatrix(8,8);
         this->ComputeElasticMatrix( HenckyElastic, ElasticMatrix);
         this->ComputePlasticHardeningParameter(HenckyElastic, Alpha, H);

         // EQUATION 1.
         AuxMatrix = Gamma * prod( AuxiliarDerivatives.PlasticPotentialDD, ElasticMatrix);
         for (unsigned int i = 0; i < 6; i++)
         {
            SystemMatrix(i,i) += 1.0;
            SystemMatrix(i,7) += AuxiliarDerivatives.PlasticPotentialD(i);
            for (unsigned int j = 0; j < 6; j++)
            {
               SystemMatrix(i,j)  += AuxMatrix(i,j);
            }
         }

         // EQUATION 2. 
         SystemMatrix(6,6) += 1.0;
         SystemMatrix(6,7) += - AuxP;

         // EQUATION 3
         AuxVector = prod( AuxiliarDerivatives.YieldFunctionD, ElasticMatrix);
         for (unsigned int i = 0; i<6; i++)
            SystemMatrix(7, i) += AuxVector(i);
         SystemMatrix(7,6) = -H;


         InverseSystemMatrix = ZeroMatrix(8,8);
         bool singular = SolidMechanicsMathUtilities<double>::InvertMatrix( SystemMatrix, InverseSystemMatrix);
         if ( singular )
         {
            std::cout << " SINGULAR MATRIX " << std::endl;
            std::cout << " SYSTEM MATRIX " << SystemMatrix << std::endl;
            std::cout << " AuxMatrix " << AuxMatrix << std::endl;
            std::cout << " ELASTIC MAT " << ElasticMatrix << std::endl;
            std::cout << " gamma " << Gamma  << std::endl;
            std::cout << " AuxiliarDerivati " << AuxiliarDerivatives.PlasticPotentialDD << std::endl;
            std::cout << " " << std::endl;
            SystemMatrix = ZeroMatrix(8,8);
            for (int i = 0; i < 8; i++)
               SystemMatrix(i,i) = 10000000.0;
            singular = SolidMechanicsMathUtilities<double>::InvertMatrix( SystemMatrix, InverseSystemMatrix);
            if ( singular )
               for (int i = 0; i < 8; i++)
                  SystemMatrix(i,i) = 0.000001;
         }

         deltaX = prod( InverseSystemMatrix, Residual);  // Compute the New Thing.

         PreviousResidualNorm = ResidualNorm;

      }

      this->CalculateKirchhoffStressVector( HenckyElastic , rStressVector);
      rNewElasticLeftCauchyGreen = ConvertHenckyStrainToCauchyGreenTensor( HenckyElastic);
      rReturnMappingVariables.DeltaGamma = Alpha; 
      rReturnMappingVariables.DeltaBeta = Gamma;  
      // CONVERGED
      double AlphaUpdate = 0.0;
      for (unsigned int j = 0; j < 3; ++j)
         AlphaUpdate += Gamma*AuxiliarDerivatives.PlasticPotentialD(j);

      double aux = 0;
      for (unsigned int j = 0; j < 3; ++j)
         aux += pow( Gamma*AuxiliarDerivatives.PlasticPotentialD(j) - AlphaUpdate/3.0, 2);

      for (unsigned int j = 3; j < 6; ++j)
         aux += 2.0*pow(Gamma*AuxiliarDerivatives.PlasticPotentialD(j) / 2.0, 2);

      rReturnMappingVariables.IncrementalPlasticShearStrain += sqrt( aux);

      /* if (iter > 60) { // STUPID FUNCION
         std::cout << "   RRMM ITER " << iter << " NORM " << ResidualNorm << std::endl;
         {
         std::cout << " CHECK THE DERIVATIVE " << std::endl;
         Vector KirS ;
         this->CalculateKirchhoffStressVector( HenckyElastic, KirS);
         Vector StressV = KirS;
         Vector NumDer = ZeroVector(6);
         double delta = 1e-5;
         double Ryield; 
         Ryield = mpYieldCriterion->CalculateYieldCondition( Ryield, StressV, Alpha);

         for (unsigned int m = 0; m < 6; m++) {
         double yield ;
         StressV = KirS;
         StressV(m) += delta; 
         yield = mpYieldCriterion->CalculateYieldCondition( yield, StressV, Alpha);
         NumDer(m) = (yield-Ryield) / delta;
         }
         std::cout << " ANALYTICAL " << AuxiliarDerivatives.YieldFunctionD << std::endl;
         std::cout << " NUMERICAL  " << NumDer << std::endl;
         std::cout << " AlphaTrial " << AlphaTrial << " ALPHA " << Alpha << std::endl;
         }
         Vector KirS ;
         this->CalculateKirchhoffStressVector( HenckyElastic, KirS);
         Vector StressV = KirS;
         Vector NumDer = ZeroVector(6);
         double delta = 1e-5;
         double Ryield; 
         Ryield = mpYieldCriterion->CalculateYieldCondition( Ryield, StressV, Alpha);

         for (unsigned int m = 0; m < 6; m++) {
         double yield ;
         StressV = KirS;
         StressV(m) += delta; 
         yield = mpYieldCriterion->CalculateYieldCondition( yield, StressV, Alpha);
         NumDer(m) = (yield-Ryield) / delta;
         }
         std::cout << " ANALYTICAL " << AuxiliarDerivatives.YieldFunctionD << std::endl;
         std::cout << " NUMERICAL  " << NumDer << std::endl;
         std::cout << " " << std::endl;
         std::cout << " " << std::endl;
      // AND NOW THE SECOND DERIVATIVE.....
      delta = 1e-5;
      for (unsigned int m = 0; m < 6; m++)
      {
      Vector Der;
      StressV = KirS;
      StressV(m) += delta; 
      mpYieldCriterion->CalculateYieldFunctionDerivative( StressV, Der, Alpha);
      for (unsigned int kk = 0; kk < 6; kk++)
      Der(kk) = (Der(kk) - AuxiliarDerivatives.YieldFunctionD(kk) ) / delta;

      std::cout << " comp " << m ;
      for ( int tt = 0; tt < 6; tt ++)
      std::cout << " " << Der(tt) ;
      std::cout << std::endl;

      }
      std::cout << " ANALYTICAL " << std::endl;;
      for (unsigned t = 0; t < 6; t++) {
      for (unsigned q = 0; q < 6; q++) {
      std::cout << " "  <<AuxiliarDerivatives.PlasticPotentialDD(t, q);
      }
      std::cout << " " << std::endl;
      }

      std::cout << " " << std::endl;
      std::cout << " " << std::endl;
      std::cout << " " << std::endl;
      std::cout << " " << std::endl;

   }*/

   }


   // THIS FUNCTION IS FOR THE RETURN MAPPING ALGORITHM (IMPLICIT)
   // IT PERFORMS SOME OPERATIONS AND GOES TO ReturnToYielSurface4
   bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMappingImpl(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
   {
      bool PlasticityActive = false;
      double Tolerance = 1e-4;

      InternalVariables PlasticVariables = mInternalVariables;
      rReturnMappingVariables.DeltaGamma = mInternalVariables.EquivalentPlasticStrain;
      rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;

      rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

      rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans( rDeltaDeformationGradient) );
      rNewElasticLeftCauchyGreen = prod( rDeltaDeformationGradient, rNewElasticLeftCauchyGreen);

      Vector StressVector;
      this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, StressVector);

      double ElasticTrialStateFunction;
      ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, StressVector, PlasticVariables.EquivalentPlasticStrain);

      if ( ElasticTrialStateFunction > Tolerance)  {
         PlasticityActive = true;
         rReturnMappingVariables.Options.Set(PLASTIC_REGION, true);
         this->ReturnStressToYieldSurface4( rReturnMappingVariables, rNewElasticLeftCauchyGreen, StressVector, ElasticTrialStateFunction, Tolerance);
      }

      rStressMatrix = MathUtils<double>::StressVectorToTensor(StressVector);

      rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
      rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

      return PlasticityActive;

   }


   // EXPLICIT STRESS INTEGRATION WITH SUBSTEPPING AND ERROR CONTROL
   // BASED ON SLOAN ET AL 2002
   bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMappingExpl(RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
   {

      //1- Initialize some variables
      bool PlasticityActive = false;
      rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

      InternalVariables PlasticVariables = mInternalVariables;

      // IncrementalPlasticStrain and Total VolumetricStrain 
      rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;
      rReturnMappingVariables.DeltaBeta = 0;
      rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain; 

      Matrix PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
      Vector PreviousStressVector;

      Vector NewStressVector;

      double Tolerance = 1.0E-4; // LMV


      //2- Check for yield Condition (at the beggining)
      this->CalculateKirchhoffStressVector( PreviousElasticLeftCauchyGreen, PreviousStressVector);
      rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, rReturnMappingVariables.DeltaGamma );

      if (rReturnMappingVariables.TrialStateFunction > Tolerance) {
         this->ReturnStressToYieldSurface( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, rReturnMappingVariables.TrialStateFunction, Tolerance);
         //THE INITIAL TENSIONAL STATE MAY BE OUTSIDE THE YIELD SURFACE?? (yes and no: no if you don't do nothing strange, however, if you apply LSInterpolation or change constant parameters during the simulation YES)
      }
      rReturnMappingVariables.DeltaBeta = 0;

      //3. Compute an Elastic Trial State
      ExplicitStressUpdateInformation  StressUpdateInformation;

      this->CalculateOneExplicitStep(rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, false, StressUpdateInformation);

      double ElasticTrialStateFunction;
      ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);

      //4a. COMPLETELY ELASTIC STEP
      if ( ElasticTrialStateFunction <= Tolerance)   {
         PlasticityActive = false;
      }

      //4b. SOME PART OF THE STEP IS PLASTIC
      else {

         PlasticityActive = true; 

         if ( ( rReturnMappingVariables.TrialStateFunction < -Tolerance) && (ElasticTrialStateFunction > Tolerance) ) {
            //2b.1 Transition from elastic to plastic
            //std::cout << " ELASTIC AND PLASTIC STEP " << std::endl;
            this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);

         }
         else {         

            bool UnloadingCondition;
            UnloadingCondition =  this->EvaluateElastoPlasticUnloadingCondition( UnloadingCondition, PreviousElasticLeftCauchyGreen, rIncrementalDeformationGradient, PlasticVariables, Tolerance); 

            if (UnloadingCondition ) {
               //2c.1 ElastoPlastic Unloading
               //std::cout << " ELASTIC UNLOADING STEP " << std::endl;
               this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);

            }
            else {
               //std::cout << " PLASTIC STEP " << std::endl;
               //2c. 2 Completedly ELastoPlastic Step 
               this->CalculateExplicitSolution( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, true, Tolerance);

            }
         }


      }


      if (PlasticityActive) {

         double DriftViolation;

         DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, rReturnMappingVariables.DeltaGamma);

         //rReturnMappingVariables.EigenVectors = NewElasticLeftCauchyGreen; 
         double Beta = rReturnMappingVariables.DeltaBeta;
         if ( fabs(DriftViolation) > Tolerance ) {

            this->ReturnStressToYieldSurface( rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, DriftViolation, Tolerance);
         }
         // LMV
         rReturnMappingVariables.DeltaBeta = Beta;


      }

      rStressMatrix = MathUtils<double>::StressVectorToTensor(NewStressVector);

      rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
      rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);


      return PlasticityActive;

   }


   bool& NonAssociativeExplicitPlasticFlowRule::EvaluateElastoPlasticUnloadingCondition( bool& rUnloadingCondition, const Matrix& rElasticLeftCauchyGreen, const Matrix& rDeltaDeformationGradient, const InternalVariables& rPlasticVariables, const double& rTolerance)
   {

      rUnloadingCondition = false;

      Vector StressVector;
      Vector YieldFunctionD;
      this->CalculateKirchhoffStressVector( rElasticLeftCauchyGreen, StressVector);
      mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, YieldFunctionD, rPlasticVariables.EquivalentPlasticStrain);

      Matrix ElasticMatrix;
      Vector ElasticStrainVector ;
      Vector IncrementalHencky;

      ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen);
      IncrementalHencky = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient) ));
      this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

      Vector IncrementalElasticStress;
      IncrementalElasticStress = prod(ElasticMatrix, IncrementalHencky);

      double Numerador = 0.0;
      double Denom1 = 0.0;
      double Denom2 = 0.0;

      for (unsigned int i = 0; i < 6; ++i) {
         Numerador += YieldFunctionD(i)*IncrementalElasticStress(i);
         Denom1 += pow( YieldFunctionD(i), 2);
         Denom2 += pow( IncrementalElasticStress(i), 2);
      }

      Denom1 = pow( Denom1*Denom2, 0.5);

      Numerador /= Denom1;

      if (Numerador < rTolerance) 
         rUnloadingCondition = true;

      return rUnloadingCondition;

   }





   // ******************   CORRECT DRIFT : CUTTING PLANE ALGORITHM ************************
   // *************************************************************************************
   void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rDrift, const double& rTolerance)
   {

      AuxiliarDerivativesStructure  AuxiliarDerivatives;
      Matrix ElasticMatrix;
      double H;
      double DeltaGamma;
      Vector ElasticCorrection;
      Matrix  CorrectedLeftCauchyGreen;
      Matrix  UpdateMatrix;
      Vector ActualElasticHenckyStrain;
      Vector StressVector;
      double Alpha = rReturnMappingVariables.DeltaGamma; 
      rReturnMappingVariables.DeltaBeta = 0.0;

      for (unsigned int i = 0; i < 150; ++i) {

         ActualElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);

         this->UpdateDerivatives(ActualElasticHenckyStrain, AuxiliarDerivatives, Alpha);
         this->ComputeElasticMatrix(ActualElasticHenckyStrain, ElasticMatrix);
         this->ComputePlasticHardeningParameter(ActualElasticHenckyStrain, Alpha, H);

         DeltaGamma = rDrift;
         DeltaGamma /= ( H + MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)));



         Vector MuyAuxiliar = -DeltaGamma*AuxiliarDerivatives.PlasticPotentialD/ 2.0;
         UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor( MuyAuxiliar);


         CorrectedLeftCauchyGreen =  prod((UpdateMatrix), rNewElasticLeftCauchyGreen);
         CorrectedLeftCauchyGreen =  prod( CorrectedLeftCauchyGreen, trans(UpdateMatrix));

         double AlphaUpdate = 0.0;
         for (unsigned int j = 0; j < 3; ++j)
            AlphaUpdate += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j);

         double aux = 0;
         for (unsigned int j = 0; j < 3; ++j)
            aux += pow(DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j) - AlphaUpdate/3.0, 2);

         for (unsigned int j = 3; j < 6; ++j)
            aux += 2.0*pow(DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j) / 2.0, 2);

         rReturnMappingVariables.IncrementalPlasticShearStrain += sqrt( aux);



         Alpha += AlphaUpdate;
         rReturnMappingVariables.DeltaBeta += DeltaGamma ;
         this->CalculateKirchhoffStressVector( CorrectedLeftCauchyGreen, rStressVector);

         rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, Alpha);

         rNewElasticLeftCauchyGreen = CorrectedLeftCauchyGreen;
         rReturnMappingVariables.DeltaGamma = Alpha;
         if( fabs(rDrift) < 1.0*rTolerance) {
            return;
         }

      }

      if ( fabs(rDrift) > 100.0*rTolerance) {
         std::cout<< " " << std::endl;
         std::cout << "Leaving drift correction WITHOUT converging " << rDrift << std::endl;
         std::cout << " StressVector " << rStressVector << std::endl;
         // COMPTUING THE LODE ANGLE TO SEE WHAT HAPPENS
       Matrix SM = ZeroMatrix(3,3);

         SM = MathUtils<double>::StressVectorToTensor( rStressVector);

         double p = 0.0;
         for (unsigned int i = 0; i < 3; ++i)
            p += SM(i,i) / 3.0;
         double J2 = 0.0;
         for (unsigned int i = 0; i < 3; ++i) {
            SM(i,i) -= p;     
            J2 += pow( SM(i,i)- p, 2);
         }
         for (unsigned int i = 0; i < 6; ++i)
            J2 += 2.0*pow( rStressVector(i), 2);
         J2 = sqrt(J2 / 2.0);

         double Lode = MathUtils<double>::Det( SM);
         Lode *= 3.0*sqrt(3.0) / 2.0;
         Lode /= pow( J2, 3);

         std::cout << " LODE " << std::asin( Lode) / 3.0 / 3.14159265359*180.0 << std::endl;
      }


   }



   //********************* UPDATES ONE ELASTOPLASTIC STRAIN INCREMENT THE STRESS ****************
   //********************************************************************************************

   void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitPlasticStep(const Matrix& rDeltaDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, const double& rPreviousEquivalentPlasticStrain, Matrix& rNewElasticLeftCauchyGreen, double& rNewEquivalentPlasticStrain, double& rNewPlasticShearStrain, double& rDeltaPlastic)
   {
      rNewEquivalentPlasticStrain = rPreviousEquivalentPlasticStrain;
      Matrix ElasticMatrix;
      Vector ElasticStrainVector ;

      ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rPreviousElasticLeftCauchyGreen);
      this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

      AuxiliarDerivativesStructure AuxiliarDerivatives;
      this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rNewEquivalentPlasticStrain);

      double H;
      this->ComputePlasticHardeningParameter(ElasticStrainVector, rNewEquivalentPlasticStrain, H);


      ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient)) );


      Vector auxVector;
      auxVector = prod(ElasticMatrix, ElasticStrainVector) ;
      double DeltaGamma = 0.0;
      for (unsigned int i = 0; i<6; ++i)
         DeltaGamma += auxVector(i)*AuxiliarDerivatives.YieldFunctionD(i);


      double auxDenominador = H + MathUtils<double>::Dot( AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD));

      DeltaGamma /= auxDenominador;

      if (DeltaGamma < 0 )
         DeltaGamma =0;

      rDeltaPlastic = DeltaGamma;

      Vector MuyAuxiliar;
      MuyAuxiliar = -DeltaGamma * AuxiliarDerivatives.PlasticPotentialD / 2.0; 

      Matrix UpdateMatrix;
      UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor (MuyAuxiliar);
      UpdateMatrix = prod( rDeltaDeformationGradient, UpdateMatrix);

      rNewElasticLeftCauchyGreen = prod(UpdateMatrix, rPreviousElasticLeftCauchyGreen);
      rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix));

      double PlasticPotentialP = 0.0;
      for (unsigned int i = 0; i < 3; ++i)
         PlasticPotentialP += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i);
      rNewEquivalentPlasticStrain +=  PlasticPotentialP ;// 3.0; 

      rNewPlasticShearStrain = 0.0;
      for (unsigned int i = 0; i < 3; ++i)
         rNewPlasticShearStrain += pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i) - PlasticPotentialP/3.0, 2);

      for (unsigned int i = 3; i < 6; ++i)
         rNewPlasticShearStrain += 2.0*pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i)/2.0, 2);

      rNewPlasticShearStrain = sqrt( rNewPlasticShearStrain );


   }

   //** EXPLICIT PART: Calculates one elasto-plastic increment with two different discretizations and computes the error  ****************
   //** it may also compute one elastic step *********************************************************************************************
   //*************************************************************************************************************************************
   void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitStep(const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, ExplicitStressUpdateInformation& rStressUpdateInformation)
   {

      if ( rElastoPlasticBool)  {

         double NewEquivalentPlasticStrain;
         double NewPlasticShearStrain;
         double DeltaPlastic; 
         Vector FirstApproxStressVector;
         // ONE STEP UPDATE
         this->CalculateOneExplicitPlasticStep( rDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables.DeltaGamma, rNewElasticLeftCauchyGreen, NewEquivalentPlasticStrain, NewPlasticShearStrain, DeltaPlastic);

         this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, FirstApproxStressVector); 

         //TWO STEP UPDATE

         Matrix HalfStepIncrementalDeformation;
         ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, 0.5, HalfStepIncrementalDeformation);

         Matrix SecondApproxLeftCauchyGreen; 

         this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables.DeltaGamma, SecondApproxLeftCauchyGreen, NewEquivalentPlasticStrain, NewPlasticShearStrain, DeltaPlastic);

         ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.5, 1.0, HalfStepIncrementalDeformation);
         rStressUpdateInformation.IncrementPlasticStrain = NewPlasticShearStrain; 
         rStressUpdateInformation.DeltaDeltaPlastic = DeltaPlastic;

         this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, SecondApproxLeftCauchyGreen, NewEquivalentPlasticStrain, rNewElasticLeftCauchyGreen, NewEquivalentPlasticStrain, NewPlasticShearStrain, DeltaPlastic);

         rStressUpdateInformation.IncrementPlasticStrain += NewPlasticShearStrain; 
         rStressUpdateInformation.DeltaDeltaPlastic += DeltaPlastic;
         rStressUpdateInformation.NewEquivalentPlasticStrain = NewEquivalentPlasticStrain;
         this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector);

         rStressUpdateInformation.StressErrorMeasure = 0.0;
         double Denominador = 0.0;
         for (unsigned int i = 0; i<6; ++i) {
            rStressUpdateInformation.StressErrorMeasure += pow( rNewStressVector(i) - FirstApproxStressVector(i), 2);
            Denominador  += pow(rNewStressVector(i), 2);
         }
         rStressUpdateInformation.StressErrorMeasure = 1.0* pow( rStressUpdateInformation.StressErrorMeasure/Denominador, 0.5);


      } 
      else  {

         rNewElasticLeftCauchyGreen = prod( rPreviousElasticLeftCauchyGreen, trans( rDeformationGradient) );
         rNewElasticLeftCauchyGreen = prod( rDeformationGradient, rNewElasticLeftCauchyGreen );

         rStressUpdateInformation.StressErrorMeasure = 0.0;

         this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector);

      }


   }






   // EXPLICIT PART: Compute the Elasto-plastic part with adaptive substepping
   void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolution( const Matrix& rIncrementalDeformationGradient, const Matrix& rPreviousElasticCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, const double& rTolerance)
   {

      KRATOS_TRY

      double TimeStep = 0.5;
      double MinTimeStep = 1.0e-4;
      double DoneTimeStep = 0.0;
      double MaxTimeStep = 0.5;
      TimeStep = MaxTimeStep;

      Matrix ActualElasticLeftCauchyGreen = rPreviousElasticCauchyGreen;
      Matrix SubstepDeformationGradient; 


      ExplicitStressUpdateInformation  StressUpdateInformation;
      bool MayBeLast = false;

      while (DoneTimeStep < 1.0)  {
         //std::cout << " TIME " << DoneTimeStep << " " << TimeStep  << std::endl;
         if (DoneTimeStep + TimeStep >= 1.0) {
            TimeStep = 1.0 - DoneTimeStep;
            MayBeLast = true;
         }

         this->ComputeSubstepIncrementalDeformationGradient( rIncrementalDeformationGradient, DoneTimeStep, DoneTimeStep + TimeStep, SubstepDeformationGradient);

         this->CalculateOneExplicitStep( SubstepDeformationGradient, ActualElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  rElastoPlasticBool, StressUpdateInformation);


         if ( StressUpdateInformation.StressErrorMeasure < rTolerance ) {
            // Se acepta el paso
            //std::cout << "   SUbstepping " << DoneTimeStep << " dt " << TimeStep << "Error " << StressErrorMeasure << std::endl;
            ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;

            this->UpdateRadialReturnVariables( rReturnMappingVariables, StressUpdateInformation);

            //double StateFunction = 0;
            /*StateFunction = mpYieldCriterion->CalculateYieldCondition( StateFunction, rNewStressVector, rReturnMappingVariables.DeltaGamma);
              if ( fabs( StateFunction) > rTolerance) {
            //this->ReturnStressToYieldSurface( rReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector, StateFunction, rTolerance);
            }*/

            DoneTimeStep += TimeStep;
            if (MayBeLast == true)
               return;

         }
         else {
            if (TimeStep <= MinTimeStep) {
               if ( StressUpdateInformation.StressErrorMeasure > rTolerance) {
                  std::cout << "ExplicitStressIntegrationDidNotConverged: StressError: " << StressUpdateInformation.StressErrorMeasure << " dt " << TimeStep <<  " meanStress " << (rNewStressVector(0) + rNewStressVector(1) + rNewStressVector(2) ) / 3.0 << std::endl ;
                  std::cout << "    stress " << rNewStressVector << std::endl;
                  double Hola;
                  Hola = mpYieldCriterion->CalculateYieldCondition(Hola, rNewStressVector, StressUpdateInformation.NewEquivalentPlasticStrain);
                  std::cout << "    rYieldValue " << Hola << std::endl;

               }

               DoneTimeStep += TimeStep;

               ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
               this->UpdateRadialReturnVariables( rReturnMappingVariables, StressUpdateInformation);
               if (MayBeLast == true )
                  return;
            }


         }

         TimeStep *= 0.9* (pow( rTolerance / (StressUpdateInformation.StressErrorMeasure+ 1e-8), 0.5 ));
         TimeStep = std::max(TimeStep, MinTimeStep);
         TimeStep = std::min(TimeStep, MaxTimeStep);
      }


      KRATOS_CATCH("")
   }





   // EXPLICIT PART:  Search for the yield surface intersection (bisection method) and compute the rest with elasto-plastic
   void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolutionWithChange(const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector,  const double& rTolerance)
   {

      // bisecction in elastic Regime

      double InitialPosition = 0.0;
      double EndPosition = 1.0;
      double HalfPosition;

      double InitialStateFunction = -1.0;
      double EndStateFunction = 1.0;
      double HalfPositionStateFunction;

      Matrix HalfPositionDeformationGradient;
      Matrix HalfPositionElasticCauchyGreen;

      ExplicitStressUpdateInformation StressUpdateInformation;

      //double NewEquivalentPlasticStrain = rReturnMappingVariables.train ;
      //double NewPlasticShearStrain  = 0.0;
      double ErrorMeasure1;
      double ErrorMeasure2;



      for (unsigned int i = 0 ; i< 150; ++i) {

         HalfPosition = 0.5*(InitialPosition + EndPosition);

         this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, HalfPosition, HalfPositionDeformationGradient);

         this->CalculateOneExplicitStep( HalfPositionDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables, HalfPositionElasticCauchyGreen, rNewStressVector, false, StressUpdateInformation);  

         HalfPositionStateFunction = mpYieldCriterion->CalculateYieldCondition( HalfPositionStateFunction, rNewStressVector, rReturnMappingVariables.DeltaGamma);


         if ( HalfPositionStateFunction < 0.0)  {
            InitialStateFunction = HalfPositionStateFunction;
            InitialPosition = HalfPosition;
         }
         else {
            EndPosition = HalfPosition;
            EndStateFunction = HalfPositionStateFunction;
         }
         ErrorMeasure1 = fabs(InitialStateFunction-EndStateFunction);
         ErrorMeasure2 = fabs(InitialPosition - EndPosition);

         if ( (ErrorMeasure1 < rTolerance ) && (ErrorMeasure2 < rTolerance)) {
            break;
         }

      }
      //std::cout << " Half Position Found " << HalfPosition << " SF " << HalfStateFunction << std::endl;

      // COMPUTE ELASTIC STEP;

      this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, HalfPosition, HalfPositionDeformationGradient);

      this->CalculateOneExplicitStep( HalfPositionDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables, HalfPositionElasticCauchyGreen, rNewStressVector, false, StressUpdateInformation);

      // COMPUTE ELASTOPLASTIC STEP

      this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, HalfPosition, 1, HalfPositionDeformationGradient); 

      this->CalculateExplicitSolution( HalfPositionDeformationGradient, HalfPositionElasticCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  true, rTolerance);


   }

   // ********************* UPDATE INTERNAL PLASTIC VARIABLES *****************************
   bool NonAssociativeExplicitPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables & rReturnMappingVariables)
   {
      mInternalVariables.EquivalentPlasticStrain = rReturnMappingVariables.DeltaGamma;

      mInternalVariables.EquivalentPlasticStrainOld = rReturnMappingVariables.IncrementalPlasticShearStrain;

      mInternalVariables.DeltaPlasticStrain += rReturnMappingVariables.IncrementalPlasticShearStrain; 

      mPlasticMultiplierVelocity = rReturnMappingVariables.DeltaBeta / rReturnMappingVariables.DeltaTime; 

      return 0;
   }

   // THIS FUNCTION IS FROM ME
   void NonAssociativeExplicitPlasticFlowRule::UpdateRadialReturnVariables( RadialReturnVariables& rReturnMappingVariables, const ExplicitStressUpdateInformation& rStressUpdateInformation)
   {
      rReturnMappingVariables.DeltaGamma = rStressUpdateInformation.NewEquivalentPlasticStrain;
      rReturnMappingVariables.IncrementalPlasticShearStrain += rStressUpdateInformation.IncrementPlasticStrain;
      rReturnMappingVariables.DeltaBeta += rStressUpdateInformation.DeltaDeltaPlastic;
   }


   void NonAssociativeExplicitPlasticFlowRule::UpdateDerivatives(const Vector& rHenckyStrain, AuxiliarDerivativesStructure& rAuxiliarDerivatives, const double& rAlpha)
   {
      Vector StressVector = ZeroVector(6);

      this->CalculateKirchhoffStressVector( rHenckyStrain, StressVector);

      this->CalculatePlasticPotentialDerivatives(StressVector, rAuxiliarDerivatives.PlasticPotentialD, rAuxiliarDerivatives.PlasticPotentialDD);

      mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, rAuxiliarDerivatives.YieldFunctionD, rAlpha);

      // in order to program things only once
      if ( rAuxiliarDerivatives.PlasticPotentialD.size() == 1 ) {
         rAuxiliarDerivatives.PlasticPotentialD = rAuxiliarDerivatives.YieldFunctionD;
         if ( rAuxiliarDerivatives.PlasticPotentialDD.size1() == 1) {
            //mpYieldCriterion->CalculateYieldFunctionSecondDerivative( StressVector, rAuxiliarDerivatives.PlasticPotentialDD, rAlpha);

            if ( rAuxiliarDerivatives.PlasticPotentialDD.size1() == 1) {
               rAuxiliarDerivatives.PlasticPotentialDD = ZeroMatrix(6,6);
            }
         }
      }
      else
      {
      } 


   }

   //***************** CalculateStress From LEFT CAUCHY GREEN MATRIX **************
   void NonAssociativeExplicitPlasticFlowRule::CalculateKirchhoffStressVector( const Matrix& rElasticLeftCauchyGreen, Vector& rNewStressVector)
   {

      Vector ElasticHenckyStrainVector;
      ElasticHenckyStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen); 
      this->CalculateKirchhoffStressVector(ElasticHenckyStrainVector, rNewStressVector);

   }

   // IDEM IN MATRIX FORM
   Matrix NonAssociativeExplicitPlasticFlowRule::ComputeKirchhoffStressMatrix( const Matrix & rElasticLeftCauchyGreen)
   {

      Vector StressVector;
      this->CalculateKirchhoffStressVector( rElasticLeftCauchyGreen, StressVector);
      Matrix StressMatrix = MathUtils<double>::StressVectorToTensor( StressVector);
      return StressMatrix;

   }

   //****** COMPUTE INCREMENTAL DEFORMATION GRADIENT  ********* 
   // FROM BATHE NOTATION ^{rFinalConfiguration} _{rReferenceConfiguration} X  from  ^{t+1}_{t} X
   //*********************************************************
   void NonAssociativeExplicitPlasticFlowRule::ComputeSubstepIncrementalDeformationGradient( const Matrix& rIncrementalDeformationGradient, const double& rReferenceConfiguration, const double& rFinalConfiguration, Matrix& rSubstepIncrementalDeformationGradient)
   {


      Matrix DeformationGradientReference;
      Matrix DeformationGradientFinal;

      Matrix IdentityMatrix = ZeroMatrix(3,3);
      for (unsigned int i = 0; i < 3; ++i)
         IdentityMatrix(i,i) = 1.0;   

      DeformationGradientReference = rReferenceConfiguration*rIncrementalDeformationGradient  + (1.0 - rReferenceConfiguration)*IdentityMatrix;
      DeformationGradientFinal     =     rFinalConfiguration*rIncrementalDeformationGradient  + (1.0 -     rFinalConfiguration)*IdentityMatrix;

      double Det; 
      MathUtils<double>::InvertMatrix(DeformationGradientReference, rSubstepIncrementalDeformationGradient, Det);
      rSubstepIncrementalDeformationGradient = prod( DeformationGradientFinal, rSubstepIncrementalDeformationGradient);

   }

   // ************** CONVERT HENCKY STRAIN TO CAUCHY-GREEN TENSOR  *********
   // **********************************************************************

   Vector NonAssociativeExplicitPlasticFlowRule::ConvertCauchyGreenTensorToHenckyStrain(const Matrix& rCauchyGreenMatrix)
   {
      Matrix EigenVectors;
      Vector EigenValues;
      int a = rCauchyGreenMatrix.size1();
      int b = rCauchyGreenMatrix.size2();
      SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues, 1e-12, 100);

      Matrix Aux = ZeroMatrix(3,3);
      a = Aux.size1();
      b = Aux.size2();
      if ( ( a != 3) | ( b != 3) )
      {
         for (int i = 0; i < 1000; i++)
         {
            std::cout << "kidding me?? already? AuxMatrix" << std::endl;
         }
      }
      for (unsigned int i = 0; i < 3; ++i)
         Aux(i,i) = (std::log(EigenValues(i)))/2.0;

      Aux = prod(Aux, (EigenVectors));
      Aux = prod(trans(EigenVectors), Aux);
      Vector Result = MathUtils<double>::StrainTensorToVector(Aux, 6);
      return Result;
   }

   // ********************** AND THE CONTRARY ***************
   Matrix NonAssociativeExplicitPlasticFlowRule::ConvertHenckyStrainToCauchyGreenTensor(const Vector& rElasticHenckyStrain)
   {

      Matrix HenckyStrainMatrix = MathUtils<double>::StrainVectorToTensor(rElasticHenckyStrain);

      Matrix EigenVectors = ZeroMatrix(3,3);
      Vector EigenValues = ZeroVector(3);

      SolidMechanicsMathUtilities<double>::EigenVectors(HenckyStrainMatrix, EigenVectors, EigenValues, 1e-12, 100);

      Matrix Aux = ZeroMatrix(3,3);

      for (unsigned int i = 0; i < 3; ++i)
         Aux(i,i) = std::exp(2.0*EigenValues(i));

      Aux = prod(Aux, (EigenVectors));
      Aux = prod(trans(EigenVectors), Aux);

      return Aux;

   }



   // ***********  ELASTO PLASTIC TANGENT MATRIX  *******
   void NonAssociativeExplicitPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rLeftCauchyGreenMatrix, const double& rAlpha, Matrix& rElasticMatrix)
   {

      int StressIntTechnique = 0; //mpYieldCriterion->GetHardeningLaw().GetProperties()[STRESS_INTEGRATION_STRATEGY];

      // COMPUTE STRAIN
      Vector ElasticStrainVector = ZeroVector(6);
      ElasticStrainVector = this->ConvertCauchyGreenTensorToHenckyStrain(rLeftCauchyGreenMatrix);

      // COMPUTE ELASTIC MATRIX
      this->ComputeElasticMatrix(ElasticStrainVector, rElasticMatrix);

      //if ( StressIntTechnique == 3) // mIMPLEX // this case allways goes with the elastic matrix.
      if ( rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE) ) {
         return; 
      }
      else if ( StressIntTechnique == 4) // IMPLEX // this case may go with elastic or the "second derivative plus"
      {
         if ( fabs( rReturnMappingVariables.NormIsochoricStress) < 1e-8 ) // ONLY ELASTIC. No Plastic part
         {
            return;
         }

         AuxiliarDerivativesStructure AuxiliarDerivatives;
         this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rAlpha);

         Matrix SystemMatrix = rElasticMatrix;

         bool singular = SolidMechanicsMathUtilities<double>::InvertMatrix( rElasticMatrix, SystemMatrix);
         if (  ! singular)
         {
            SystemMatrix += rReturnMappingVariables.NormIsochoricStress *  AuxiliarDerivatives.PlasticPotentialDD;

            singular = SolidMechanicsMathUtilities<double>::InvertMatrix( SystemMatrix, SystemMatrix);
         }
         if ( !singular)
            rElasticMatrix = SystemMatrix;

         return;
      }

      // NOW WE ARE IN THE IMPLICIT OR EXPLICIT.

      if ( rReturnMappingVariables.Options.Is( FlowRule::PLASTIC_REGION ) )
      {

         AuxiliarDerivativesStructure AuxiliarDerivatives;
         this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rAlpha);


         // ADD the second derivative in the Elastic Matrix
         if ( StressIntTechnique == 2) // IMPLICIT CASE
         {
            Matrix SystemMatrix = rElasticMatrix;

            bool singular = SolidMechanicsMathUtilities<double>::InvertMatrix( rElasticMatrix, SystemMatrix);
            if ( ! singular) {
               SystemMatrix += rReturnMappingVariables.DeltaBeta *  AuxiliarDerivatives.PlasticPotentialDD;

               singular = SolidMechanicsMathUtilities<double>::InvertMatrix( SystemMatrix, SystemMatrix);
            }
            if ( !singular)
            {
               rElasticMatrix = SystemMatrix;
            }
         }

         // MATRIX 1.

         double H;
         this->ComputePlasticHardeningParameter(ElasticStrainVector, rAlpha, H);

         Vector AuxVectorF;
         Vector AuxVectorG;

         AuxVectorF = prod(trans(AuxiliarDerivatives.YieldFunctionD), rElasticMatrix);
         AuxVectorG = prod( rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD);

         Matrix PlasticUpdate;
         PlasticUpdate = MyCrossProduct(rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD, AuxiliarDerivatives.YieldFunctionD);

         rElasticMatrix -= 1.0*PlasticUpdate / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));

      }


   }


   // IS LINE SEARCH... however // vale, esto no sirve para nada, lo programé porque no entendia una cosa del implex, pero ya lo vi claro
   void NonAssociativeExplicitPlasticFlowRule::PerformSomeSortOfLineSearch( Vector & rHenckyElastic, double& rGamma, double& rAlpha, const Vector & rHenckyTrialElastic, const double& rAlphaTrial, double & ResidualNorm, const double& rPreviousError, const Vector & rDeltaX)
   {

      Vector HenckyStrain = rHenckyElastic;
      double Gamma = rGamma; 
      double Alpha = rAlpha;

      // desandar lo andado
      for (int i = 0; i < 6; i++)
         HenckyStrain(i) += rDeltaX(i);

      Alpha += rDeltaX(6);
      Gamma += rDeltaX(7);

      double beta, BestBeta = 0;
      double Ratio = 0;
      AuxiliarDerivativesStructure  AuxiliarDerivatives;

      Vector Hen, Res = ZeroVector(8), R1 = ZeroVector(6), StressV;
      double R2, R3, AuxP; 
      double Alp, Gam;


      for ( int arc = -50; arc < 100; arc++)
      {
         if ( arc == 0)
            continue;

         beta = double( arc) / 50.0;

         Hen = HenckyStrain;
         Alp = Alpha; 
         Gam = Gamma; 

         // ACTUALIZAR

         for ( int i = 0; i < 6; i++)
            Hen(i) -= beta*rDeltaX(i);

         Alp -= beta * rDeltaX(6);
         Gam -= beta * rDeltaX(7);


         // AND NOW I HAVE TO COMPUTE THE ERROR...
         this->UpdateDerivatives( Hen, AuxiliarDerivatives, Alp);

         R1 = Hen + Gam * AuxiliarDerivatives.PlasticPotentialD - rHenckyTrialElastic;

         AuxP = 0;
         for ( int i = 0; i < 3; i++)
            AuxP += AuxiliarDerivatives.PlasticPotentialD(i);

         R2 = Alp - rAlphaTrial - Gam * AuxP;

         this->CalculateKirchhoffStressVector( Hen, StressV);
         R3 = mpYieldCriterion->CalculateYieldCondition( R3, StressV, Alp);

         for ( int i = 0; i < 6; i++)
            Res(i) = R1(i);
         Res(6) = R2; 
         Res(7) = R3;

         ResidualNorm = 0;
         for ( int i = 0; i < 8; i++)
            ResidualNorm += pow ( Res(i) , 2);
         ResidualNorm = sqrt( ResidualNorm);

         if ( rPreviousError / ResidualNorm  > Ratio)
         {
            Ratio = rPreviousError / ResidualNorm;
            BestBeta = beta;
         }

         //std::cout << " THis Line Search: " << beta << " bestB " << BestBeta << " RATIO " << rPreviousError / ResidualNorm << " bestRatio " << Ratio << std::endl;
         //std::cout << rPreviousError << " and " << ResidualNorm << std::endl;
      }

      ResidualNorm =  rPreviousError / Ratio;
      for (int i = 0; i < 6; i++)
         rHenckyElastic(i) = HenckyStrain(i) - BestBeta * rDeltaX(i);

      rAlpha = Alpha - BestBeta * rDeltaX(6);
      rGamma = Gamma - BestBeta * rDeltaX(7);
      //std::cout << " BETA BEST " << BestBeta << std::endl;
      //if (Ratio < 1)
      //std::cout << " this Line did not converge " << std::endl;

   }

   Matrix NonAssociativeExplicitPlasticFlowRule::MyCrossProduct(const Matrix& rM, const Vector& rA, const Vector& rB)
   {


      Vector A = rA;
      Vector B = rB;

      A = prod(rM, A);
      B = prod(trans(B), rM); 

      Matrix Result = ZeroMatrix(6,6);
      for (unsigned int i = 0; i<6; ++i) {
         for ( unsigned int j = 0; j<6; ++j) {
            Result(i,j) = A(i)*B(j);
         }
      }

      return Result;

   }

   //  *********   SAVE FUNCTIONS ***********
   void NonAssociativeExplicitPlasticFlowRule::save( Serializer& rSerializer ) const
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
         rSerializer.save("mPlasticMultiplierVelocity",mPlasticMultiplierVelocity);
   }

   void NonAssociativeExplicitPlasticFlowRule::load( Serializer& rSerializer )
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
         rSerializer.load("mPlasticMultiplierVelocity",mPlasticMultiplierVelocity);
   }



} //End Namepace Kratos
