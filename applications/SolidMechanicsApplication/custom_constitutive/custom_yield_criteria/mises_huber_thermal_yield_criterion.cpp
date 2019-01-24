//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_constitutive/custom_yield_criteria/mises_huber_thermal_yield_criterion.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

MisesHuberThermalYieldCriterion::MisesHuberThermalYieldCriterion()
	:MisesHuberYieldCriterion()
{

}


//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

MisesHuberThermalYieldCriterion::MisesHuberThermalYieldCriterion(HardeningLawPointer pHardeningLaw)
	:MisesHuberYieldCriterion(pHardeningLaw)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MisesHuberThermalYieldCriterion& MisesHuberThermalYieldCriterion::operator=(MisesHuberThermalYieldCriterion const& rOther)
{
   MisesHuberYieldCriterion::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

MisesHuberThermalYieldCriterion::MisesHuberThermalYieldCriterion(MisesHuberThermalYieldCriterion const& rOther)
	:MisesHuberYieldCriterion(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

YieldCriterion::Pointer MisesHuberThermalYieldCriterion::Clone() const
{
  return Kratos::make_shared<MisesHuberThermalYieldCriterion>(*this);
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

MisesHuberThermalYieldCriterion::~MisesHuberThermalYieldCriterion()
{
}

/// Operations.


//***************************CALCULATE PLASTIC DISSIPATION****************************
//************************************************************************************

double& MisesHuberThermalYieldCriterion::CalculatePlasticDissipation(double & rPlasticDissipation, const Parameters& rValues)
{
        const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const HardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double Hardening = 0;

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

double& MisesHuberThermalYieldCriterion::CalculateDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues)
{

	const double& rLameMu_bar = rValues.GetLameMu_bar();
	const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const HardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double DeltaHardening = 0;

	// std::cout<<" delta dissipation A "<<std::endl;
	// rHardeningParameters.print();

	DeltaHardening = mpHardeningLaw->CalculateDeltaHardening( DeltaHardening, rHardeningParameters );

	double Hardening = 0;


	// std::cout<<" delta dissipation B "<<std::endl;
	Hardening = mpHardeningLaw->CalculateHardening( Hardening, rHardeningParameters );

	double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

       	double DeltaThermalHardening = 0;

	// std::cout<<" delta dissipation C "<<std::endl;
	//rHardeningParameters.print();

	DeltaThermalHardening = mpHardeningLaw->CalculateDeltaThermalHardening( DeltaThermalHardening, rHardeningParameters );

	// std::cout<<"  DeltaHardening  "<<DeltaHardening<<" DeltaThermalHardening "<<DeltaThermalHardening<<std::endl;
	// std::cout<<"  EquivalentStress  "<<EquivalentStress<<" DeltaGamma "<<rDeltaGamma<<std::endl;

	rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);

	rDeltaPlasticDissipation *= ( (-1) * DeltaThermalHardening );

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

double& MisesHuberThermalYieldCriterion::CalculateImplexPlasticDissipation(double & rPlasticDissipation, const Parameters& rValues)
{

	const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const HardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();

	double Hardening = 0;

	Hardening = mpHardeningLaw->CalculateHardening( Hardening, rHardeningParameters );

	//PENDENT(change the definition of this stress Hardening has a different expression  !!!!)
	double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

	rPlasticDissipation = 0.9 * EquivalentStress * rDeltaGamma * ( 1.0/rDeltaTime );

	return rPlasticDissipation;
}


//*****************CALCULATE IMPLEX DELTA PLASTIC DISSIPATION*************************
//************************************************************************************

double& MisesHuberThermalYieldCriterion::CalculateImplexDeltaPlasticDissipation(double & rDeltaPlasticDissipation, const Parameters& rValues)
{

	double DeltaThermalHardening = 0;

	const double& rDeltaGamma = rValues.GetDeltaGamma();
	const double& rDeltaTime  = rValues.GetDeltaTime();

	const HardeningLaw::Parameters& rHardeningParameters = rValues.GetHardeningParameters();


	DeltaThermalHardening = mpHardeningLaw->CalculateDeltaThermalHardening( DeltaThermalHardening, rHardeningParameters );


	rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);

	rDeltaPlasticDissipation *= ( (-1) * DeltaThermalHardening );

	rDeltaPlasticDissipation *= rDeltaGamma;

	return rDeltaPlasticDissipation;
}


void MisesHuberThermalYieldCriterion::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MisesHuberYieldCriterion );
}

void MisesHuberThermalYieldCriterion::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MisesHuberYieldCriterion );
}


}  // namespace Kratos.
