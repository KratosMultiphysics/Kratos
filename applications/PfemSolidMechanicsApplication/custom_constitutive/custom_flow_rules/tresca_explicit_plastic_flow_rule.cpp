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
#include "custom_constitutive/custom_flow_rules/tresca_explicit_plastic_flow_rule.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{



//************ CONSTRUCTOR ***********
TrescaExplicitFlowRule::TrescaExplicitFlowRule()
   :J2ExplicitFlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

TrescaExplicitFlowRule::TrescaExplicitFlowRule(YieldCriterionPointer pYieldCriterion)
	:J2ExplicitFlowRule(pYieldCriterion)
{
   
}

//********* ASSIGNMENT OPERATOR
TrescaExplicitFlowRule& TrescaExplicitFlowRule::operator=(TrescaExplicitFlowRule const& rOther)
{
	J2ExplicitFlowRule::operator=(rOther);
	return *this;

}



//********** COPY CONSTRUCTOR *********
TrescaExplicitFlowRule::TrescaExplicitFlowRule(TrescaExplicitFlowRule const& rOther)
      :J2ExplicitFlowRule(rOther)
{
}

//*******   CLONE ********
FlowRule::Pointer TrescaExplicitFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new TrescaExplicitFlowRule(*this));
  return p_clone;
}



// ********** DESTRUCTOR **************
TrescaExplicitFlowRule::~TrescaExplicitFlowRule()
{
}





void TrescaExplicitFlowRule::ComputePlasticHardeningParameter(const Vector& rHenckyStrainVector, const double& rAlpha, double& rH)
{
   
     rH = 0.0;

}
 
void TrescaExplicitFlowRule::CalculatePlasticPotentialDerivatives(const Vector& rStressVector, Vector& rFirstDerivative, Matrix& rSecondDerivative)
{
     //double YieldStress = mpYieldCriterion->GetHardeningLaw().GetProperties()[YIELD_STRESS];
     rFirstDerivative = ZeroVector(1);
     rSecondDerivative = ZeroMatrix(1,1);
     return;

}
void TrescaExplicitFlowRule::save( Serializer& rSerializer) const 
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, J2ExplicitFlowRule )
}

void TrescaExplicitFlowRule::load( Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, J2ExplicitFlowRule )

}

} //end namespace kratos
