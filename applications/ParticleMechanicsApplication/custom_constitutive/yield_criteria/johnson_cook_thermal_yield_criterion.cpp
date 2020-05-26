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

/// Operations.


//***************************CALCULATE YIELD CONDITION********************************
//************************************************************************************

double& JohnsonCookThermalYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Parameters& rValues)
{
	double Hardening = 0.0;

	const double& rStressNorm = rValues.GetStressNorm();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	//std::cout<<" yield function "<<std::endl;
	//rHardeningParameters.print();

	Hardening = mpHardeningLaw->CalculateHardening(Hardening, rHardeningParameters);

	rStateFunction = rStressNorm - sqrt(2.0 / 3.0) * Hardening;

	return rStateFunction;
}


//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& JohnsonCookThermalYieldCriterion::CalculateStateFunction(double& rStateFunction, const Parameters& rValues)
{

	const double& rStressNorm = rValues.GetStressNorm();
	const double& rLameMu_bar = rValues.GetLameMu_bar();
	const double& rDeltaGamma = rValues.GetDeltaGamma();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();


	double Hardening = 0.0;

	//std::cout<<" state function "<<std::endl;
	//rHardeningParameters.print();

	Hardening = mpHardeningLaw->CalculateHardening(Hardening, rHardeningParameters);

	rStateFunction = rStressNorm - 2.0 * rLameMu_bar * rDeltaGamma - sqrt(2.0 / 3.0) * (Hardening);

	return rStateFunction;
}


//***************************CALCULATE STATE FUNCTION ********************************
//************************************************************************************

double& JohnsonCookThermalYieldCriterion::CalculateDeltaStateFunction(double& rDeltaStateFunction, const Parameters& rValues)
{
	const double& rLameMu_bar = rValues.GetLameMu_bar();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double DeltaHardening = 0.0;


	//std::cout<<" delta state function "<<std::endl;
	//rHardeningParameters.print();

	DeltaHardening = mpHardeningLaw->CalculateDeltaHardening(DeltaHardening, rHardeningParameters);

	//std::cout<<" DeltaHardening "<<DeltaHardening<<std::endl;

	rDeltaStateFunction = 2.0 * rLameMu_bar + (2.0 / 3.0) * DeltaHardening;

	return rDeltaStateFunction;
}


//***************************CALCULATE PLASTIC DISSIPATION****************************
//************************************************************************************

double& JohnsonCookThermalYieldCriterion::CalculatePlasticDissipation(double & rPlasticDissipation, const Parameters& rValues)
{
        const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double Hardening = 0.0;

	//std::cout<<" dissipation "<<std::endl;
	//rHardeningParameters.print();

	Hardening = mpHardeningLaw->CalculateHardening( Hardening, rHardeningParameters );


	double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

	rPlasticDissipation = 0.9 * EquivalentStress * rDeltaGamma * ( 1.0/rDeltaTime );

	//std::cout<<"  Hardening  "<<Hardening<<" rDeltaGamma "<<rDeltaGamma<<" EquivalentStress "<<EquivalentStress<<" PlasticDissipation "<<rPlasticDissipation<<std::endl;


	return rPlasticDissipation;
}


//**********************CALCULATE DELTA PLASTIC DISSIPATION***************************
//************************************************************************************

double& JohnsonCookThermalYieldCriterion::CalculateDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues)
{

	const double& rLameMu_bar = rValues.GetLameMu_bar();
	const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double DeltaHardening = 0.0;

	// std::cout<<" delta dissipation A "<<std::endl;
	// rHardeningParameters.print();

	DeltaHardening = mpHardeningLaw->CalculateDeltaHardening( DeltaHardening, rHardeningParameters );

	double Hardening = 0.0;


	// std::cout<<" delta dissipation B "<<std::endl;
	Hardening = mpHardeningLaw->CalculateHardening( Hardening, rHardeningParameters );

	double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

       	double DeltaThermalHardening = 0.0;

	// std::cout<<" delta dissipation C "<<std::endl;
	//rHardeningParameters.print();

	DeltaThermalHardening = mpHardeningLaw->CalculateDeltaThermalHardening( DeltaThermalHardening, rHardeningParameters );

	// std::cout<<"  DeltaHardening  "<<DeltaHardening<<" DeltaThermalHardening "<<DeltaThermalHardening<<std::endl;
	// std::cout<<"  EquivalentStress  "<<EquivalentStress<<" DeltaGamma "<<rDeltaGamma<<std::endl;

	rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);

	rDeltaPlasticDissipation *= ( (-1.0) * DeltaThermalHardening );

	rDeltaPlasticDissipation *= (rDeltaGamma - ( EquivalentStress + DeltaHardening * rDeltaGamma * (2.0/3.0) )/( 2.0 * rLameMu_bar + (2.0/3.0) * DeltaHardening ) );

	//std::cout<<" Temperature "<<Temperature<<" ReferenceTemperature "<<ReferenceTemperature<<std::endl;
	//std::cout<<" DeltaThermalHardening "<<DeltaThermalHardening<<std::endl;
	//std::cout<<" DeltaGamma "<<rDeltaGamma<<std::endl;
	//std::cout<<" DeltaHardening "<<DeltaHardening<<std::endl;

	// //Revise Dissipation in Johnson-Cook	material

	// //----Dependent plastic rate
	// rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);
	// rDeltaPlasticDissipation *= (-1) * rDeltaThermalHardening;
	// rDeltaPlasticDissipation *= ( rDeltaGamma - (EquivalentStress + (2.0/3.0) * (rDeltaGamma) * DeltaHardening)/ (2.0 * LameMu_bar + (2.0/3.0) * DeltaHardening);

	// //----Independent plastic rate
	// rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);
	// rDeltaPlasticDissipation *= (-1) * rDeltaThermalHardening;
	// rDeltaPlasticDissipation *= ( rDeltaGamma - (Equivalentstress + (2.0/3.0) * rDeltaHardening * rDeltaGamma)/(2.0 * LameMu_bar + (2.0/3.0) * DeltaHardening + Coef6);

	// Coef6 = ( (2.0/3.0) * ( A + B * pow ( rVariables.EquivalentPlasticStrain, n) ) * NormalizedTemperature * C / (sqrt(2.0/3.0) * rParameters.DeltaGamma) ) ;

        // //Revise Dissipation in Baker Johnson-Cook	material
	// //----Dependent plastic rate  and //----Independent plastic rate
	// rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);
	// rDeltaPlasticDissipation *= (-1) * rDeltaThermalHardening;

	// rDeltaPlasticDissipation *= ( rDeltaGamma - (EquivalentStress + (2.0/3.0) * (rDeltaGamma) * DeltaHardening)/ (2.0 * LameMu_bar + (2.0/3.0) * DeltaHardening);


	return rDeltaPlasticDissipation;
}




//***********************CALCULATE IMPLEX PLASTIC DISSIPATION*************************
//************************************************************************************

double& JohnsonCookThermalYieldCriterion::CalculateImplexPlasticDissipation(double & rPlasticDissipation, const Parameters& rValues)
{

	const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double Hardening = 0;

	Hardening = mpHardeningLaw->CalculateHardening( Hardening, rHardeningParameters );

	//PENDENT(change the definition of this stress Hardening has a different expression  !!!!)
	double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

	rPlasticDissipation = 0.9 * EquivalentStress * rDeltaGamma * ( 1.0/rDeltaTime );

	return rPlasticDissipation;
}


//*****************CALCULATE IMPLEX DELTA PLASTIC DISSIPATION*************************
//************************************************************************************

double& JohnsonCookThermalYieldCriterion::CalculateImplexDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues)
{

	double DeltaThermalHardening = 0.0;

	const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const ParticleHardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();


	DeltaThermalHardening = mpHardeningLaw->CalculateDeltaThermalHardening( DeltaThermalHardening, rHardeningParameters );


	rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);

	rDeltaPlasticDissipation *= ( (-1.0) * DeltaThermalHardening );

	rDeltaPlasticDissipation *= rDeltaGamma;

	return rDeltaPlasticDissipation;
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
