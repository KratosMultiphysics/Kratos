//
//   Project Name:        KratosPoromechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/continuum_laws/custom_yield_criteria/modified_mises_yield_criterion.hpp"

#include "poromechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

ModifiedMisesYieldCriterion::ModifiedMisesYieldCriterion()
	:YieldCriterion()
{

}


//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

ModifiedMisesYieldCriterion::ModifiedMisesYieldCriterion(HardeningLawPointer pHardeningLaw)
	:YieldCriterion(pHardeningLaw)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ModifiedMisesYieldCriterion& ModifiedMisesYieldCriterion::operator=(ModifiedMisesYieldCriterion const& rOther)
{
   YieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

ModifiedMisesYieldCriterion::ModifiedMisesYieldCriterion(ModifiedMisesYieldCriterion const& rOther)
	:YieldCriterion(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

YieldCriterion::Pointer ModifiedMisesYieldCriterion::Clone() const
{
  return Kratos::make_shared<ModifiedMisesYieldCriterion>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ModifiedMisesYieldCriterion::~ModifiedMisesYieldCriterion()
{
}

/// Operations.


//************************** CALCULATE EQUIVALENT STRAIN *****************************
//************************************************************************************

double& ModifiedMisesYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Parameters& rVariables)
{
    // Compute I1
    const Matrix& StrainMatrix = rVariables.GetStrainMatrix();
    const unsigned int Dim = StrainMatrix.size1();
    double I1 = 0.0;

    for(unsigned int i = 0; i < Dim; i++)
    {
        I1 += StrainMatrix(i,i);
    }

    // Compute J2
    Matrix DeviatoricStrain(Dim,Dim);
    noalias(DeviatoricStrain) = StrainMatrix;

    for(unsigned int i = 0; i < Dim; i++)
    {
        DeviatoricStrain(i,i) -= I1/Dim;
    }

    Matrix Auxiliar(Dim,Dim);
    noalias(Auxiliar) = prod(DeviatoricStrain,DeviatoricStrain);
    double J2 = 0.0;

    for(unsigned int i = 0; i < Dim; i++)
    {
        J2 += Auxiliar(i,i);
    }

    J2 *= 0.5;

    // Compute Equivalent Strain (rStateFunction)
    const Properties& MaterialProperties = mpHardeningLaw->GetProperties();
    const double& StrengthRatio = MaterialProperties[STRENGTH_RATIO];
    const double& PoissonRatio = MaterialProperties[POISSON_RATIO];

    rStateFunction = I1*(StrengthRatio-1.0)/(2.0*StrengthRatio*(1.0-2.0*PoissonRatio)) +
                    sqrt( I1*I1*(StrengthRatio-1.0)*(StrengthRatio-1.0)/((1.0-2.0*PoissonRatio)*(1.0-2.0*PoissonRatio)) +
                          J2*12.0*StrengthRatio/((1.0+PoissonRatio)*(1.0+PoissonRatio)) )/(2.0*StrengthRatio);

    return rStateFunction;
}


//***************************CALCULATE DAMAGE PARAMETER ******************************
//************************************************************************************

double& ModifiedMisesYieldCriterion::CalculateStateFunction(double& rStateFunction, const Parameters& rVariables)
{
    const HardeningLaw::Parameters& HardeningLawParameters = rVariables.GetHardeningParameters();

    mpHardeningLaw->CalculateHardening(rStateFunction, HardeningLawParameters);

	return rStateFunction;
}


//***************************CALCULATE DAMAGE DERIVATIVE *****************************
//************************************************************************************

double& ModifiedMisesYieldCriterion::CalculateDeltaStateFunction(double& rDeltaStateFunction, const Parameters& rVariables)
{
    const HardeningLaw::Parameters& HardeningLawParameters = rVariables.GetHardeningParameters();

    mpHardeningLaw->CalculateDeltaHardening(rDeltaStateFunction, HardeningLawParameters);

	return rDeltaStateFunction;
}


//************************************************************************************
//************************************************************************************

void ModifiedMisesYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void ModifiedMisesYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
