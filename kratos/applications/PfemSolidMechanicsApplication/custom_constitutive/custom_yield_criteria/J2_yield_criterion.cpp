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

   double YieldStress = this->GetHardeningLaw().GetProperties()[YIELD_STRESS];
 
   rStateFunction = 0.0;

   double MeanStress = 0.0;
  
   for (unsigned int i = 0; i<3; ++i)
       MeanStress += rStressVector(i);

   MeanStress /= 3.0;

   for (unsigned int i = 0; i<3; ++i)
       rStateFunction += pow(rStressVector(i) - MeanStress, 2.0);

   for (unsigned int i = 3; i<6; ++i)
       rStateFunction += 2.0*pow(rStressVector(i), 2.0);


   rStateFunction = pow(3.0/2.0, 1.0/2.0)*pow( rStateFunction , 1.0/2.0);

   // ASSUMING THAT YIELD STRESS IS EQUAL TO S_u

   rStateFunction -= 2.0*YieldStress;;
//   rStateFunction = -1.0;
   return rStateFunction; 
}




//************************* YIELD FUNCTION DERIVATIVE ******************
//**********************************************************************
void J2YieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rYieldFunctionD, const double& rAlpha)
{


     double MeanStress = 0.0;
 
     for (unsigned int i = 0; i<3; ++i)
         MeanStress += rStressVector(i);

     MeanStress /= 3.0;

     double denominador = 0.0;
     Vector ShearVector = ZeroVector(6);

     for (unsigned int i = 0; i<3; ++i)  {
         ShearVector(i) = rStressVector(i) - MeanStress; 
         denominador += pow( ShearVector(i), 2.0);
     }

    for (unsigned int i = 3; i<6; ++i) {
         ShearVector(i) = 2.0*rStressVector(i);
         denominador += 2.0*pow( ShearVector(i), 2.0);
     }

     denominador = pow(3.0/2.0, 1.0/2.0)*pow(denominador, 1.0/2.0);

     rYieldFunctionD = ShearVector;
 
     rYieldFunctionD *= 3.0/2.0/denominador; 



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
