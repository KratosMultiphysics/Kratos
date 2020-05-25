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

// External includes

// Project includes
#include "custom_constitutive/yield_criteria/mises_huber_yield_criterion.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

MisesHuberYieldCriterion::MisesHuberYieldCriterion()
	:ParticleYieldCriterion()
{

}


//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MisesHuberYieldCriterion::MisesHuberYieldCriterion(HardeningLawPointer pHardeningLaw)
	:ParticleYieldCriterion(pHardeningLaw)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MisesHuberYieldCriterion& MisesHuberYieldCriterion::operator=(MisesHuberYieldCriterion const& rOther)
{
	ParticleYieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

MisesHuberYieldCriterion::MisesHuberYieldCriterion(MisesHuberYieldCriterion const& rOther)
	:ParticleYieldCriterion(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ParticleYieldCriterion::Pointer MisesHuberYieldCriterion::Clone() const
{
  return Kratos::make_shared<MisesHuberYieldCriterion>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

MisesHuberYieldCriterion::~MisesHuberYieldCriterion()
{
}

/// Operations.


//***************************CALCULATE YIELD CONDITION********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateYieldCondition(double & rStateFunction, const Parameters& rValues)
{
	double Hardening = 0.0;

	const double& rStressNorm = rValues.GetStressNorm();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	//std::cout<<" yield function "<<std::endl;
	//rHardeningParameters.print();

	Hardening = mpHardeningLaw->CalculateHardening(Hardening, rHardeningParameters);

	rStateFunction = rStressNorm - sqrt(2.0/3.0) * Hardening;

	return rStateFunction;
}


//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateStateFunction(double & rStateFunction, const Parameters& rValues)
{

	const double& rStressNorm = rValues.GetStressNorm();
	const double& rLameMu_bar = rValues.GetLameMu_bar();
	const double& rDeltaGamma = rValues.GetDeltaGamma();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();


	double Hardening = 0.0;

	//std::cout<<" state function "<<std::endl;
	//rHardeningParameters.print();

	Hardening = mpHardeningLaw->CalculateHardening( Hardening, rHardeningParameters );

	rStateFunction = rStressNorm - 2.0 * rLameMu_bar * rDeltaGamma - sqrt(2.0/3.0) * ( Hardening );

	return rStateFunction;
}


//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateDeltaStateFunction(double & rDeltaStateFunction, const Parameters& rValues)
{
	const double& rLameMu_bar = rValues.GetLameMu_bar();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double DeltaHardening = 0.0;


	//std::cout<<" delta state function "<<std::endl;
	//rHardeningParameters.print();

	DeltaHardening = mpHardeningLaw->CalculateDeltaHardening( DeltaHardening, rHardeningParameters );

	//std::cout<<" DeltaHardening "<<DeltaHardening<<std::endl;

	rDeltaStateFunction = 2.0 * rLameMu_bar + (2.0/3.0) * DeltaHardening;

	return rDeltaStateFunction;
}


//***************************CALCULATE PLASTIC DISSIPATION****************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculatePlasticDissipation(double & rPlasticDissipation, const Parameters& rValues)
{
	rPlasticDissipation = 0.0;
	return rPlasticDissipation;
}


//**********************CALCULATE DELTA PLASTIC DISSIPATION***************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues)
{
	rDeltaPlasticDissipation = 0.0;
	return rDeltaPlasticDissipation;
}


//***************************CALCULATE PLASTIC DISSIPATION****************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateImplexPlasticDissipation(double & rPlasticDissipation, const Parameters& rValues)
{
	rPlasticDissipation = 0.0;
	return rPlasticDissipation;
}


//**********************CALCULATE DELTA PLASTIC DISSIPATION***************************
//************************************************************************************

double& MisesHuberYieldCriterion::CalculateImplexDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues)
{
	rDeltaPlasticDissipation = 0.0;
	return rDeltaPlasticDissipation;
}



//************************************************************************************
//************************************************************************************

void MisesHuberYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ParticleYieldCriterion )
}

void MisesHuberYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ParticleYieldCriterion )
}


}  // namespace Kratos.
