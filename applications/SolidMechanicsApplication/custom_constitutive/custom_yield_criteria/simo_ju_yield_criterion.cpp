//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solid_mechanics_application_variables.h"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

SimoJuYieldCriterion::SimoJuYieldCriterion()
	:YieldCriterion()
{
   
}


//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

SimoJuYieldCriterion::SimoJuYieldCriterion(HardeningLawPointer pHardeningLaw)
	:YieldCriterion(pHardeningLaw)
{
   
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SimoJuYieldCriterion& SimoJuYieldCriterion::operator=(SimoJuYieldCriterion const& rOther)
{
   YieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

SimoJuYieldCriterion::SimoJuYieldCriterion(SimoJuYieldCriterion const& rOther)
	:YieldCriterion(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

YieldCriterion::Pointer SimoJuYieldCriterion::Clone() const
{
  YieldCriterion::Pointer p_clone(new SimoJuYieldCriterion(*this));
  return p_clone;
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SimoJuYieldCriterion::~SimoJuYieldCriterion()
{
}

/// Operations.


//************************** CALCULATE EQUIVALENT STRAIN *****************************
//************************************************************************************

double& SimoJuYieldCriterion::CalculateYieldCondition(double & rStateFunction, const Parameters& rVariables)
{
    const Properties& MaterialProperties = mpHardeningLaw->GetProperties();
    const double& StrengthRatio = MaterialProperties[STRENGTH_RATIO];
    
    const Matrix& StrainMatrix = rVariables.GetStrainMatrix();
    const Matrix& StressMatrix = rVariables.GetStressMatrix();

    Vector PrincipalStresses = SolidMechanicsMathUtilities<double>::EigenValues(StressMatrix,1e-10,1e-10);

    double Macaulay_PrincipalStress,Absolute_PrincipalStress;
    double Theta = 0.0,Theta_den = 0.0;

    for(unsigned int i=0;i<PrincipalStresses.size();i++)
    {
        Macaulay_PrincipalStress = PrincipalStresses[i];
        Absolute_PrincipalStress = PrincipalStresses[i];
        if(Macaulay_PrincipalStress<0.0) Macaulay_PrincipalStress = 0.0;
        if(Absolute_PrincipalStress<0.0) Absolute_PrincipalStress = -Absolute_PrincipalStress;
        Theta += Macaulay_PrincipalStress;
        Theta_den += Absolute_PrincipalStress;
    }

    if(Theta_den > 1e-20) Theta = Theta/Theta_den;
    else Theta = 0.5;

    Matrix Auxiliar = prod(StrainMatrix,StressMatrix);
    double StressNorm = 0.0;
    for(unsigned int i=0; i<Auxiliar.size1(); i++)
        StressNorm += Auxiliar(i,i);

    rStateFunction = (Theta+(1.0-Theta)/StrengthRatio)*sqrt(StressNorm);
    
    return rStateFunction;
}


//***************************CALCULATE DAMAGE PARAMETER ******************************
//************************************************************************************

double& SimoJuYieldCriterion::CalculateStateFunction(double& rStateFunction, const Parameters& rVariables)
{    
    const HardeningLaw::Parameters& HardeningLawParameters = rVariables.GetHardeningParameters();
    
    mpHardeningLaw->CalculateHardening(rStateFunction, HardeningLawParameters);
    
	return rStateFunction;
}


//***************************CALCULATE DAMAGE DERIVATIVE *****************************
//************************************************************************************

double& SimoJuYieldCriterion::CalculateDeltaStateFunction(double& rDeltaStateFunction, const Parameters& rVariables)
{
    const HardeningLaw::Parameters& HardeningLawParameters = rVariables.GetHardeningParameters();
    
    mpHardeningLaw->CalculateDeltaHardening(rDeltaStateFunction, HardeningLawParameters);
    
	return rDeltaStateFunction;
}


//************************************************************************************
//************************************************************************************

void SimoJuYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void SimoJuYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
