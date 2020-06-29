//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

// System includes

// External includes

// Project includes

#include "custom_constitutive/rht_concrete_3D_law.hpp"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
	RHTConcrete3DLaw::RHTConcrete3DLaw()
		: HyperElastic3DLaw()
	{  }


	RHTConcrete3DLaw::RHTConcrete3DLaw(const RHTConcrete3DLaw& rOther)
		: HyperElastic3DLaw(rOther)
	{  }


	ConstitutiveLaw::Pointer RHTConcrete3DLaw::Clone() const
	{
		return Kratos::make_shared<RHTConcrete3DLaw>(*this);
	}


	RHTConcrete3DLaw::~RHTConcrete3DLaw()
	{  }


	bool RHTConcrete3DLaw::Has(const Variable<double>& rThisVariable)
	{
		if (rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN
			|| rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN_RATE
			|| rThisVariable == MP_EQUIVALENT_STRESS)
			return true;
		else return false;
	}


	void RHTConcrete3DLaw::InitializeMaterial(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
	{
		BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
	}


	void RHTConcrete3DLaw::CalculateMaterialResponseKirchhoff(Kratos::ConstitutiveLaw::Parameters& rValues)
	{
		// Check if the constitutive parameters are passed correctly to the law calculation
		CheckParameters(rValues);

		// Get Values to compute the constitutive law:
		Flags& Options = rValues.GetOptions();
		KRATOS_ERROR_IF(Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
			<< "The JohnsonCookThermalPlastic3DLaw cannot accept the deformation graddient F as a strain input."
			<< " Please set the strain vector and set USE_ELEMENT_PROVIDED_STRAIN = False in the element.";
		const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
		CheckIsExplicitTimeIntegration(CurrentProcessInfo);
		const Properties& MaterialProperties = rValues.GetMaterialProperties();
		const double yield_stress_failure_ratio = 1e-3; // particle fails if current_yield/virgin_yield < yield_stress_failure_ratio

		// Get old stress vector and current strain vector
		const Vector StrainVector = rValues.GetStrainVector();
		Vector& StressVector = rValues.GetStressVector();
		const Matrix identity = (GetStrainSize() == 3) ? IdentityMatrix(2) : IdentityMatrix(3);

		// Convert vectors to matrices for easier manipulation
		Matrix stress_old = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		Matrix strain_increment = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		//MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
		//MakeStrainStressMatrixFromVector(StressVector, stress_old);


	}


	int RHTConcrete3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

		KRATOS_ERROR_IF (JC_PARAMETER_A.Key()==0 || rMaterialProperties[JC_PARAMETER_A] < 0.0) << "JC_PARAMETER_A has key zero or invalid value (expected positive number ~500MPa)" << std::endl;
		KRATOS_ERROR_IF (JC_PARAMETER_B.Key()==0 || rMaterialProperties[JC_PARAMETER_B] < 0.0) << "JC_PARAMETER_B has key zero or invalid value (expected positive number ~500MPa)" << std::endl;
		KRATOS_ERROR_IF (JC_PARAMETER_C.Key()==0 || rMaterialProperties[JC_PARAMETER_C] < 0.0) << "JC_PARAMETER_C has key zero or invalid value (expected positive number ~0.01)" << std::endl;
		KRATOS_ERROR_IF (JC_PARAMETER_n.Key()==0 || rMaterialProperties[JC_PARAMETER_n] < 0.0) << "JC_PARAMETER_n has key zero or invalid value (expected positive number ~0.25)" << std::endl;
		KRATOS_ERROR_IF (REFERENCE_STRAIN_RATE.Key()==0 || rMaterialProperties[REFERENCE_STRAIN_RATE] <= 0.0) << "REFERENCE_STRAIN_RATE has key zero or invalid value (expected positive number ~1.0)" << std::endl;
		KRATOS_ERROR_IF (TAYLOR_QUINNEY_COEFFICIENT.Key()==0 || rMaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] < 0.0) << "TAYLOR_QUINNEY_COEFFICIENT has key zero or invalid value (expected positive number ~0.9)" << std::endl;

		if (rMaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] > 0.0)
		{
			// Check parameters that affect thermal softening
			KRATOS_ERROR_IF(JC_PARAMETER_m.Key() == 0 || rMaterialProperties[JC_PARAMETER_m] < 0.0) << "JC_PARAMETER_m has key zero or invalid value (expected positive number ~1.0)" << std::endl;
			KRATOS_ERROR_IF(MELD_TEMPERATURE.Key() == 0 || rMaterialProperties[MELD_TEMPERATURE] <= 0.0) << "MELD_TEMPERATURE has key zero or invalid value (expected positive number ~1700K)" << std::endl;
			KRATOS_ERROR_IF(REFERENCE_TEMPERATURE.Key() == 0 || rMaterialProperties[REFERENCE_TEMPERATURE] <= 0.0) << "REFERENCE_TEMPERATURE has key zero or invalid value (expected positive number ~293K)" << std::endl;
			KRATOS_ERROR_IF(TEMPERATURE.Key() == 0 || rMaterialProperties[TEMPERATURE] <= 0.0) << "TEMPERATURE has key zero or invalid value (expected positive number ~293K)" << std::endl;
			KRATOS_ERROR_IF(SPECIFIC_HEAT.Key() == 0 || rMaterialProperties[SPECIFIC_HEAT] < 0.0) << "SPECIFIC_HEAT has key zero or invalid value (expected positive number ~450.0)" << std::endl;
		}

		if (check_base > 1) return 1;
		return 0;
	}


	bool RHTConcrete3DLaw::CheckParameters(Parameters& rValues)
	{
		return rValues.CheckAllParameters();
	}


	void RHTConcrete3DLaw::MakeStrainStressVectorFromMatrix(const Matrix& rInput, Vector& rOutput)
	{
		if (rOutput.size() != GetStrainSize()) rOutput.resize(GetStrainSize(), false);

		// 3D stress arrangement
		// Normal components
		rOutput[0] = rInput(0, 0);
		rOutput[1] = rInput(1, 1);
		rOutput[2] = rInput(2, 2);

		// Shear components
		rOutput[3] = 2.0 * rInput(0, 1); //xy
		rOutput[4] = 2.0 * rInput(1, 2); //yz
		rOutput[5] = 2.0 * rInput(0, 2); //xz
	}


	void RHTConcrete3DLaw::MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput)
	{
		if (rOutput.size1() != 3 || rOutput.size2() != 3)rOutput.resize(3, 3, false);

		// 3D stress arrangement
		// Normal components
		rOutput(0, 0) = rInput[0];
		rOutput(1, 1) = rInput[1];
		rOutput(2, 2) = rInput[2];

		// Shear components
		rOutput(0, 1) = 0.5 * rInput[3]; //xy
		rOutput(1, 2) = 0.5 * rInput[4]; //yz
		rOutput(0, 2) = 0.5 * rInput[5]; //xz

		// Fill symmetry
		rOutput(1, 0) = rOutput(0, 1);
		rOutput(2, 1) = rOutput(1, 2);
		rOutput(2, 0) = rOutput(0, 2);
	}


	void RHTConcrete3DLaw::CheckIsExplicitTimeIntegration(const ProcessInfo& rCurrentProcessInfo)
	{
		const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
			? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
			: false;
		KRATOS_ERROR_IF_NOT(is_explicit) << "The Johnson Cook MPM material law is currently limited to explicit time integration only";
	}


	double RHTConcrete3DLaw::CalculateHardenedYieldStress(const Properties& MaterialProperties,
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		// Hardening formula is: = (A + B* ep^n) * strain_rate_hardening_factor * thermal_hardening_factor
		double hardened_stress = MaterialProperties[JC_PARAMETER_A] + MaterialProperties[JC_PARAMETER_B] *
			std::pow(EquivalentPlasticStrain, MaterialProperties[JC_PARAMETER_n]);
		//hardened_stress *= CalculateStrainRateHardeningFactor(MaterialProperties, PlasticStrainRate);
		//hardened_stress *= CalculateThermalHardeningFactor(MaterialProperties, Temperature);

		return hardened_stress;
	}



	void RHTConcrete3DLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
		rFeatures.mOptions.Set(FINITE_STRAINS);
		rFeatures.mOptions.Set(ISOTROPIC);

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_Velocity_Gradient);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the spacedimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}


	double& RHTConcrete3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == MP_EQUIVALENT_STRESS)
		{
			rValue = mEquivalentStress;
		}
		else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function GetValue double.";

		return(rValue);
	}


	void RHTConcrete3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function SetValue double.";
	}
} // Namespace Kratos