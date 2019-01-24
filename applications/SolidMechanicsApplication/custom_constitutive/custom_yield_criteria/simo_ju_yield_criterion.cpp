//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"

#include "solid_mechanics_application_variables.h"

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
  return Kratos::make_shared<SimoJuYieldCriterion>(*this);
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
    // Compute Theta parameter
    double Theta;

    const Matrix& StressMatrix = rVariables.GetStressMatrix();
    const unsigned int Dim = StressMatrix.size1();

    Vector PrincipalStresses(Dim);
    if(Dim == 2)
    {
        PrincipalStresses[0] = 0.5*(StressMatrix(0,0)+StressMatrix(1,1)) +
                               sqrt(0.25*(StressMatrix(0,0)-StressMatrix(1,1))*(StressMatrix(0,0)-StressMatrix(1,1)) +
                                    StressMatrix(0,1)*StressMatrix(0,1));
        PrincipalStresses[1] = 0.5*(StressMatrix(0,0)+StressMatrix(1,1)) -
                               sqrt(0.25*(StressMatrix(0,0)-StressMatrix(1,1))*(StressMatrix(0,0)-StressMatrix(1,1)) +
                                    StressMatrix(0,1)*StressMatrix(0,1));
    }
    else
    {
        noalias(PrincipalStresses) = SolidMechanicsMathUtilities<double>::EigenValuesDirectMethod(StressMatrix);
    }

    double Macaulay_PrincipalStress = 0.0, Absolute_PrincipalStress = 0.0;

    for(unsigned int i = 0; i < Dim; i++)
    {
        if(PrincipalStresses[i] > 0.0)
        {
            Macaulay_PrincipalStress += PrincipalStresses[i];
            Absolute_PrincipalStress += PrincipalStresses[i];
        }
        else
        {
            Absolute_PrincipalStress -= PrincipalStresses[i];
        }
    }

    if(Absolute_PrincipalStress > 1.0e-20)
    {
        Theta = Macaulay_PrincipalStress/Absolute_PrincipalStress;
    }
    else
    {
        Theta = 0.5;
    }

    // Compute Equivalent Strain (rStateFunction)
    const Matrix& StrainMatrix = rVariables.GetStrainMatrix();
    Matrix Auxiliar(Dim,Dim);
    noalias(Auxiliar) = prod(StrainMatrix,StressMatrix);

    double StressNorm = 0.0;

    for(unsigned int i = 0; i < Dim; i++)
    {
        StressNorm += Auxiliar(i,i);
    }

    const double& StrengthRatio = mpHardeningLaw->GetProperties()[STRENGTH_RATIO];

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
