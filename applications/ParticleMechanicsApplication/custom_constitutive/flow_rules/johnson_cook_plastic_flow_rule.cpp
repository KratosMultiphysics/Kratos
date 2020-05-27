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
#include "includes/properties.h"
#include "custom_constitutive/flow_rules/johnson_cook_plastic_flow_rule.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************

	JohnsonCookPlasticFlowRule::JohnsonCookPlasticFlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

	JohnsonCookPlasticFlowRule::JohnsonCookPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:ParticleFlowRule(pYieldCriterion)
{

}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

	JohnsonCookPlasticFlowRule& JohnsonCookPlasticFlowRule::operator=(JohnsonCookPlasticFlowRule const& rOther)
{
   ParticleFlowRule::operator=(rOther);
   return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

	JohnsonCookPlasticFlowRule::JohnsonCookPlasticFlowRule(JohnsonCookPlasticFlowRule const& rOther)
	:ParticleFlowRule(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ParticleFlowRule::Pointer JohnsonCookPlasticFlowRule::Clone() const
{
  return Kratos::make_shared<JohnsonCookPlasticFlowRule>(*this);
}


//********************************DESTRUCTOR******************************************
//************************************************************************************

JohnsonCookPlasticFlowRule::~JohnsonCookPlasticFlowRule()
{
}

/// Operations.

//***************************CALCULATE STRESS NORM ***********************************
//************************************************************************************

double& JohnsonCookPlasticFlowRule::CalculateStressNorm(Matrix& rStressMatrix, double& rStressNorm)
{
	rStressNorm = std::sqrt((rStressMatrix(0, 0) * rStressMatrix(0, 0)) + (rStressMatrix(1, 1) * rStressMatrix(1, 1)) + (rStressMatrix(2, 2) * rStressMatrix(2, 2)) +
	(rStressMatrix(0, 1) * rStressMatrix(0, 1)) + (rStressMatrix(0, 2) * rStressMatrix(0, 2)) + (rStressMatrix(1, 2) * rStressMatrix(1, 2)) +
	(rStressMatrix(1, 0) * rStressMatrix(1, 0)) + (rStressMatrix(2, 0) * rStressMatrix(2, 0)) + (rStressMatrix(2, 1) * rStressMatrix(2, 1)));


	return rStressNorm;
}

double JohnsonCookPlasticFlowRule::CalculateThermalDerivative(const ParticleHardeningLaw::Parameters& rValues)
{
	//get values
	const double rPlasticStrainRate = rValues.GetPlasticStrainRate();
	const double rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double rTemperature = rValues.GetTemperature();

	//Constant Parameters of the -- Johnson and Cook --:
	const double A = GetProperties()[JC_PARAMETER_A];
	const double B = GetProperties()[JC_PARAMETER_B];
	const double C = GetProperties()[JC_PARAMETER_C];

	const double n = GetProperties()[JC_PARAMETER_n];
	const double m = GetProperties()[JC_PARAMETER_m];

	const double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	const double MeldTemperature = GetProperties()[MELD_TEMPERATURE];
	const double ReferenceStrainRate = GetProperties()[REFERENCE_STRAIN_RATE];

	double thermal_derivative = 0.0;
	if (ReferenceTemperature <= rTemperature && rTemperature <= MeldTemperature)
	{
		double strain_rate_hardening_factor = 1.0;
		if (rPlasticStrainRate > ReferenceStrainRate)
		{
			strain_rate_hardening_factor += C * std::log(rPlasticStrainRate / ReferenceStrainRate);
		}

		double thermal_hardening_factor = std::pow((rTemperature - ReferenceTemperature) / (MeldTemperature - ReferenceTemperature), m);

		double temp = -1.0 * m * (A + B * std::pow(rEquivalentPlasticStrain, n)) / (rTemperature - ReferenceTemperature);

		thermal_derivative = temp * strain_rate_hardening_factor * thermal_hardening_factor;
	}
	return thermal_derivative;
}

double JohnsonCookPlasticFlowRule::CalculatePlasticStrainRateDerivative(const ParticleHardeningLaw::Parameters& rValues)
{
	//get values
	const double rPlasticStrainRate = rValues.GetPlasticStrainRate();
	const double rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double rTemperature = rValues.GetTemperature();

	//Constant Parameters of the -- Johnson and Cook --:
	const double A = GetProperties()[JC_PARAMETER_A];
	const double B = GetProperties()[JC_PARAMETER_B];
	const double C = GetProperties()[JC_PARAMETER_C];

	const double n = GetProperties()[JC_PARAMETER_n];
	const double m = GetProperties()[JC_PARAMETER_m];

	const double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	const double MeldTemperature = GetProperties()[MELD_TEMPERATURE];
	const double ReferenceStrainRate = GetProperties()[REFERENCE_STRAIN_RATE];

	double plastic_strain_rate_derivative = 0.0;

	if (rPlasticStrainRate >= ReferenceStrainRate)
	{
		double thermal_hardening_factor;
		if (rTemperature < ReferenceTemperature)
		{
			thermal_hardening_factor = 1.0;
		}
		else if (rTemperature >= MeldTemperature)
		{
			thermal_hardening_factor = 0.0;
		}
		else
		{
			thermal_hardening_factor = 1.0 - std::pow((rTemperature - ReferenceTemperature) / (MeldTemperature - ReferenceTemperature), m);
		}

		plastic_strain_rate_derivative = C / rPlasticStrainRate * (A + B * std::pow(rEquivalentPlasticStrain, n)) * thermal_hardening_factor;
	}

	return plastic_strain_rate_derivative;
}

double JohnsonCookPlasticFlowRule::CalculatePlasticStrainDerivative(const ParticleHardeningLaw::Parameters& rValues)
{
	//get values
	const double rPlasticStrainRate = rValues.GetPlasticStrainRate();
	const double rEquivalentPlasticStrain = rValues.GetEquivalentPlasticStrain();
	const double rTemperature = rValues.GetTemperature();

	//Constant Parameters of the -- Johnson and Cook --:
	const double A = GetProperties()[JC_PARAMETER_A];
	const double B = GetProperties()[JC_PARAMETER_B];
	const double C = GetProperties()[JC_PARAMETER_C];

	const double n = GetProperties()[JC_PARAMETER_n];
	const double m = GetProperties()[JC_PARAMETER_m];

	const double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
	const double MeldTemperature = GetProperties()[MELD_TEMPERATURE];
	const double ReferenceStrainRate = GetProperties()[REFERENCE_STRAIN_RATE];

	double plastic_strain_derivative = 0.0;

	// Calculate thermal hardening factor
	double thermal_hardening_factor;
	if (rTemperature < ReferenceTemperature)
	{
		thermal_hardening_factor = 1.0;
	}
	else if (rTemperature >= MeldTemperature)
	{
		thermal_hardening_factor = 0.0;
	}
	else
	{
		thermal_hardening_factor = 1.0 - std::pow((rTemperature - ReferenceTemperature) / (MeldTemperature - ReferenceTemperature), m);
	}

	// Calculate strain rate hardening factor
	double strain_rate_hardening_factor = 1.0;
	if (rPlasticStrainRate > ReferenceStrainRate)
	{
		strain_rate_hardening_factor += C * std::log(rPlasticStrainRate / ReferenceStrainRate);
	}

	plastic_strain_derivative = n * B * std::pow(rEquivalentPlasticStrain, (n - 1.0)) * strain_rate_hardening_factor * thermal_hardening_factor;

	return plastic_strain_derivative;
}

//***************************CALCULATE RADIAL RETURN MAPPING**************************
//************************************************************************************

bool JohnsonCookPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix)
{

	//0.- Initialize Variables
	bool PlasticityActive = false;
	rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

	InternalVariables PlasticVariables = mInternalVariables;
	ParticleYieldCriterion::Parameters CriterionParameters;
	this->SetCriterionParameters(rReturnMappingVariables, PlasticVariables, CriterionParameters);


	//1.- Isochoric stress norm
	rReturnMappingVariables.NormIsochoricStress = CalculateStressNorm(rIsoStressMatrix, rReturnMappingVariables.NormIsochoricStress);

	//2.- Check yield condition
	rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, CriterionParameters);


	//3.- Initialize PlasticDissipation
	mThermalVariables.PlasticDissipation = 0;
	mThermalVariables.DeltaPlasticDissipation = 0;



	if (rReturnMappingVariables.Options.Is(IMPLEX_ACTIVE))
	{

		this->CalculateImplexReturnMapping(rReturnMappingVariables, PlasticVariables, CriterionParameters, rIsoStressMatrix);

	}
	else {

		if (rReturnMappingVariables.TrialStateFunction <= 0)
		{

			PlasticityActive = false;
			PlasticVariables.DeltaPlasticStrain = 0;
			rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

		}
		else
		{

			//3.- Calculate the consistency condition
			bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, PlasticVariables, CriterionParameters);

			if (!converged)
				std::cout << " ConstitutiveLaw did not converge " << std::endl;


			//4.- Update back stress, plastic strain and stress
			UpdateConfiguration(rReturnMappingVariables, rIsoStressMatrix);


			//5.- Calculate thermal dissipation and delta thermal dissipation
			this->CalculateThermalDissipation(CriterionParameters, rReturnMappingVariables.Thermal);

			PlasticityActive = true;
			rReturnMappingVariables.Options.Set(PLASTIC_REGION, true);
		}

	}

	// std::cout<<" ReturnMapping "<<std::endl;
	// mInternalVariables.print();
	// mThermalVariables.print();
	// std::cout<<" rIsoStressMatrix "<<rIsoStressMatrix<<std::endl;

	rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

	return 	PlasticityActive;
}


void JohnsonCookPlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ParticleFlowRule );
}

void JohnsonCookPlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ParticleFlowRule );
}


}  // namespace Kratos.
