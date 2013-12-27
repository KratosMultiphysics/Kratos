// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_flow_rules/J2_explicit_plastic_flow_rule.hpp"
#include "utilities/math_utils.h"
#include "includes/ublas_interface.h"
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
      double Young = 2.069E+05;
      double Nu = 0.29;
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


} //end namespace kratos
