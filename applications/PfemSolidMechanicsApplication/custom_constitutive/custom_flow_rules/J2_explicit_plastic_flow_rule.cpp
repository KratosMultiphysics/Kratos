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
#include "custom_constitutive/custom_flow_rules/J2_explicit_plastic_flow_rule.hpp"


#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{



   //************ CONSTRUCTOR ***********
   J2ExplicitFlowRule::J2ExplicitFlowRule()
      :NonAssociativeExplicitPlasticFlowRule()
   {
   }

   //*****************************INITIALIZATION CONSTRUCTOR*****************************
   //************************************************************************************

   J2ExplicitFlowRule::J2ExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
      :NonAssociativeExplicitPlasticFlowRule(pYieldCriterion)
   {

   }

   //********* ASSIGMENT OPERATOR
   J2ExplicitFlowRule& J2ExplicitFlowRule::operator=(J2ExplicitFlowRule const& rOther)
   {
      NonAssociativeExplicitPlasticFlowRule::operator=(rOther);
      return *this;
   }


   //********** COPY CONSTRUCTOR *********
   J2ExplicitFlowRule::J2ExplicitFlowRule(J2ExplicitFlowRule const& rOther)
      :NonAssociativeExplicitPlasticFlowRule(rOther)
   {
   }


   //*******   CLONE ********
   FlowRule::Pointer J2ExplicitFlowRule::Clone() const
   {
      FlowRule::Pointer p_clone(new J2ExplicitFlowRule(*this));
      return p_clone;
   }


   // ********** DESTRUCTOR **************
   J2ExplicitFlowRule::~J2ExplicitFlowRule()
   {
   }



   void J2ExplicitFlowRule::CalculateKirchhoffStressVector(const Vector& rHenckyStrainVector, Vector& rKirchhoffStressVector)
   {
      Matrix ElasticMatrix;
      this->ComputeElasticMatrix(rHenckyStrainVector, ElasticMatrix);

      rKirchhoffStressVector = prod(ElasticMatrix, rHenckyStrainVector);

   }


   void J2ExplicitFlowRule::ComputeElasticMatrix(const Vector& rElasticStrainVector, Matrix& rElasticMatrix)
   {

      Matrix Aux = ZeroMatrix(6,6);
      rElasticMatrix = Aux;

      double Young = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
      double Nu = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
      double diagonal =   Young/(1.0+Nu)/(1.0-2.0*Nu) * (1.0-Nu);
      double nodiagonal = Young/(1.0+Nu)/(1.0-2.0*Nu) * ( Nu);
      double corte      = Young/(1.0+Nu)/2.0;

      for (unsigned int i = 0; i<3; ++i) {
         for (unsigned int j = 0; j<3; ++j) {
            if (i == j) {
               rElasticMatrix(i,i) = diagonal;
            }
            else {
               rElasticMatrix(i,j) = nodiagonal;
            }
         }
      }

      for (unsigned int j = 3; j<6; ++j)
         rElasticMatrix(j,j) = corte;

   }


   void J2ExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
   {

      rH = 0.0;

   }

   void J2ExplicitFlowRule::CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative) 
   {
      rFirstDerivative = ZeroVector(1);

      rSecondDerivative = ZeroMatrix(6,6);
      Vector DevStress = rStressVector; 

      double MeanStress = 0.0;
      for (unsigned int i = 0; i < 3; i++)
      {
         MeanStress += rStressVector(i) / 3.0;
      }

      for (unsigned int i = 0; i < 3; i++)
      {
         DevStress(i) -= MeanStress; 
      }

      double denominador = 0.0;
      for (unsigned int i = 0; i < 3; i++)
      {
            denominador += pow( DevStress(i) , 2) ;
      }
      for (unsigned int i = 3; i < 6; i++)
      {
            denominador += 2.0*pow( DevStress(i) , 2) ;
      }

      denominador = sqrt(denominador / 2.0);
      if ( denominador < 1e-5)
         return; 

      
      for (unsigned int i = 0; i < 6; i++) {
         for ( unsigned int j = 0; j < 6; j++) {
            double times = 0.5;
            if ( i > 2)
               times *= 2.0;
            if (j > 2)
               times *= 2.0;
            rSecondDerivative(i,j) -=  times * DevStress(i) * DevStress(j) / denominador/ denominador;
         }
      }

      for (unsigned int i = 0; i <3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            rSecondDerivative(i,j) -=  1.0/3.0;
         }
      }

      for (unsigned int i = 0 ; i < 6; i++) {
         double times = 1;
         if ( i > 2)
            times = 2; 
         rSecondDerivative(i,i) += times;
      }


      rSecondDerivative *= sqrt(3) / ( 2.0 * denominador) ;

   }
   void J2ExplicitFlowRule::save( Serializer& rSerializer) const 
   {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
   }

   void J2ExplicitFlowRule::load( Serializer& rSerializer)
   {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )

   }

} //end namespace kratos
