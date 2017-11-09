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
#include "custom_constitutive/custom_yield_criteria/J2_yield_criterion.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


//*******************************CONSTRUCTOR******************************************
//************************************************************************************
J2YieldCriterion::J2YieldCriterion()
	:YieldCriterion()
{
   
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

J2YieldCriterion::J2YieldCriterion(HardeningLawPointer pHardeningLaw)
	:YieldCriterion(pHardeningLaw)
{
   
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

J2YieldCriterion& J2YieldCriterion::operator=(J2YieldCriterion const& rOther)
{
   YieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

J2YieldCriterion::J2YieldCriterion(J2YieldCriterion const& rOther)
	:YieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

J2YieldCriterion::~J2YieldCriterion()
{
}



//************************* CALCULATE YIELD FUNCTION  ******************
//**********************************************************************

double& J2YieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha)
{
   // es la J2 = sqrt( 1/2 s_ij s_ij) multiplicat per sqrt(3);

   rStateFunction = 0.0;
   double YieldStress = this->GetHardeningLaw().GetProperties()[YIELD_STRESS];

   double MeanStress, J2;

   StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2);

   rStateFunction = sqrt(3.0) * J2 - 2.0 * YieldStress; 


   return rStateFunction; 
}




//************************* YIELD FUNCTION DERIVATIVE ******************
//**********************************************************************
void J2YieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
{

     double MeanStress, J2;
     Vector V1, V2; 

     StressInvariantsUtilities::CalculateStressInvariants( rStressVector, MeanStress, J2);
     StressInvariantsUtilities::CalculateDerivativeVectors( rStressVector, V1, V2);

     rYieldFunctionD = sqrt(3.0)*V2;

}
void J2YieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void J2YieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
