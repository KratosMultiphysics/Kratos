//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_constitutive/yield_criteria/johnson_cook_thermal_yield_criterion.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

JohnsonCookThermalYieldCriterion::JohnsonCookThermalYieldCriterion()
	:ParticleYieldCriterion()
{

}


//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

JohnsonCookThermalYieldCriterion::JohnsonCookThermalYieldCriterion(HardeningLawPointer pHardeningLaw)
	:ParticleYieldCriterion(pHardeningLaw)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

JohnsonCookThermalYieldCriterion& JohnsonCookThermalYieldCriterion::operator=(JohnsonCookThermalYieldCriterion const& rOther)
{
   ParticleYieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

JohnsonCookThermalYieldCriterion::JohnsonCookThermalYieldCriterion(JohnsonCookThermalYieldCriterion const& rOther)
	:ParticleYieldCriterion(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ParticleYieldCriterion::Pointer JohnsonCookThermalYieldCriterion::Clone() const
{
  return Kratos::make_shared<JohnsonCookThermalYieldCriterion>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

JohnsonCookThermalYieldCriterion::~JohnsonCookThermalYieldCriterion()
{
}

double& JohnsonCookThermalYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Parameters& rValues)
{
	const double& rStressNorm = rValues.GetStressNorm();

	double Hardening = 0.0;

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	Hardening = mpHardeningLaw->CalculateHardening(Hardening, rHardeningParameters);

	rStateFunction = rStressNorm - Hardening;

	return rStateFunction;
}

void JohnsonCookThermalYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ParticleYieldCriterion );
}

void JohnsonCookThermalYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ParticleYieldCriterion );
}


}  // namespace Kratos.
