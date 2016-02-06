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

      Matrix Aux = ZeroMatrix(6);
      rElasticMatrix = Aux;

      double Young      = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
      double Nu         = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
      double diagonal   = Young/(1.0+Nu)/(1.0-2.0*Nu) * (1.0-Nu);
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
    rSecondDerivative = ZeroMatrix(1);
    return;  // program things only once
     double YieldStress = mpYieldCriterion->GetHardeningLaw().GetProperties()[YIELD_STRESS]
;


     rSecondDerivative = ZeroMatrix(6);

     double MeanStress  = 0.0;
     for (unsigned int i = 0; i < 3; ++i)
         MeanStress += rStressVector(i);
     MeanStress /= 3.0;

     double Denominador = 0.0;
     Vector ShearVector = ZeroVector(6);

     for (unsigned int i = 0; i < 3; ++i) {
         Denominador += pow( rStressVector(i) - MeanStress, 2.0);
         ShearVector(i) = rStressVector(i) - MeanStress;
     }

     for (unsigned int i = 3; i < 6; ++i) {
         Denominador += 2.0 * pow ( rStressVector(i) , 2.0);
         ShearVector(i) = 2.0 * rStressVector(i);
     }

     if (Denominador < 1e-8) {
         return;
     }
     for (unsigned int i = 0; i < 3; ++i)
         rSecondDerivative(i,i) = 1.0; 

     for (unsigned int i = 3; i < 6; ++i)
         rSecondDerivative(i,i) = 2.0; 

     for (unsigned int i = 0; i < 3; ++i) {
         for (unsigned int j = 0; j < 3; ++j) {
             rSecondDerivative(i,j) = -1.0/3.0;
         }
      }

     rSecondDerivative *= sqrt (3.0/2.0) / sqrt(Denominador);

     for (unsigned int i = 0; i < 6; ++i) {
         for (unsigned int j = 0; j < 6; ++j) {
             rSecondDerivative(i,j) += - sqrt(3.0/2.0) * ShearVector(i) * ShearVector(j)  / pow( Denominador, 3.0/2.0); 
         }
     }
    
     rSecondDerivative /= YieldStress ; 

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
