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

#include "custom_constitutive/johnson_cook_thermal_plastic_3D_law.hpp"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"
#include "mpm_application_variables.h"

namespace Kratos
{
	JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw()
		: HyperElastic3DLaw()
	{  }


	JohnsonCookThermalPlastic3DLaw::JohnsonCookThermalPlastic3DLaw(const JohnsonCookThermalPlastic3DLaw& rOther)
		: HyperElastic3DLaw(rOther)
	{  }


	ConstitutiveLaw::Pointer JohnsonCookThermalPlastic3DLaw::Clone() const
	{
		return Kratos::make_shared<JohnsonCookThermalPlastic3DLaw>(*this);
	}


	JohnsonCookThermalPlastic3DLaw::~JohnsonCookThermalPlastic3DLaw()
	{  }


	bool JohnsonCookThermalPlastic3DLaw::Has(const Variable<double>& rThisVariable)
	{
		if (rThisVariable == MP_TEMPERATURE
			|| rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN
			|| rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN_RATE
			|| rThisVariable == MP_HARDENING_RATIO
			|| rThisVariable == MP_EQUIVALENT_STRESS)
			return true;
		else return false;
	}


	void JohnsonCookThermalPlastic3DLaw::InitializeMaterial(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
	{
		BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

		mStrainOld = ZeroVector(GetStrainSize());
		mEquivalentPlasticStrainOld = 0.0;
		mPlasticStrainRateOld = 0.0;
		mEnergyInternal = 0.0;
		mEnergyDissipated = 0.0;
		mTemperatureOld = rMaterialProperties[TEMPERATURE];
		mGammaOld = 1e-8;
		mHardeningRatio = 1.0;

		if (rMaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] == 0.0) {
			KRATOS_WARNING("Johnson Cook Material Model") << " Taylor Quinney Coefficient set to 0, ignoring thermal effects" << std::endl;
		}

		mYieldStressOld = CalculateHardenedYieldStress(rMaterialProperties, mEquivalentPlasticStrainOld, mPlasticStrainRateOld, mTemperatureOld);
		mYieldStressVirgin = mYieldStressOld;
	}


	void JohnsonCookThermalPlastic3DLaw::CalculateMaterialResponseKirchhoff(Kratos::ConstitutiveLaw::Parameters& rValues)
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
		MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
		MakeStrainStressMatrixFromVector(StressVector, stress_old);

		// Material moduli
		const double shear_modulus_G = MaterialProperties[YOUNG_MODULUS] / (2.0 + 2.0 * MaterialProperties[POISSON_RATIO]);
		const double bulk_modulus_K = MaterialProperties[YOUNG_MODULUS] / (3.0 - 6.0 * MaterialProperties[POISSON_RATIO]);

		// Calculate deviatoric quantities
		const Matrix strain_increment_deviatoric = strain_increment -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(strain_increment) / 3.0 * identity;
		const Matrix stress_deviatoric_old = stress_old -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_old) / 3.0 * identity;

		// Calculate trial (predicted) j2 stress
		double stress_hydrostatic_new = MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_old) / 3.0 +
			bulk_modulus_K * MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(strain_increment);
		Matrix stress_deviatoric_trial = stress_deviatoric_old + 2.0 * shear_modulus_G * strain_increment_deviatoric;
		const double j2_stress_trial = std::sqrt(3.0 / 2.0 *
			MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_deviatoric_trial));

		// Declare deviatoric stress matrix to be used later
		Matrix stress_deviatoric_converged = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);

		// Assume current yield stress is the same as the old (elastic predictor)
		double yield_stress = mYieldStressOld;

		if (j2_stress_trial > yield_stress && yield_stress/mYieldStressVirgin > yield_stress_failure_ratio)
		{
			// Newton raphson setup
			double gamma = mGammaOld;
			double gamma_min = 0.0;
			double gamma_max = j2_stress_trial / GetSqrt6() / shear_modulus_G;
			bool is_converged = false;
			const SizeType iteration_limit = 100;
			IndexType iteration = 0;
			const double tolerance_delta_gamma = 1e-9;
			double yield_function, yield_function_gradient, dYield_dGamma,
				delta_gamma, predicted_eps, predicted_eps_rate,
				predicted_temperature;

			// Flow direction
			Matrix flow_direction_normalized = stress_deviatoric_trial / (j2_stress_trial / std::sqrt(3.0 / 2.0));

			// Initial prediction of quantities
			predicted_eps = mEquivalentPlasticStrainOld + GetSqrt23() * gamma; // eps = equivalent plastic strain
			predicted_eps_rate = GetSqrt23() * gamma / CurrentProcessInfo[DELTA_TIME];
			predicted_temperature = mTemperatureOld;
			predicted_temperature += MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] / GetSqrt6() / MaterialProperties[DENSITY]
				/ MaterialProperties[SPECIFIC_HEAT] * (yield_stress + mYieldStressOld) * gamma;

			// Newton Raphson return mapping loop
			while (!is_converged)
			{
				// Calculate predicted yield stress
				yield_stress = CalculateHardenedYieldStress(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);

				// Compute yield function and derivative
				yield_function = j2_stress_trial - GetSqrt6() * shear_modulus_G * gamma - yield_stress;
				if (yield_function < 0.0) gamma_max = gamma;
				else gamma_min = gamma;
				dYield_dGamma = CalculatePlasticStrainDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
				dYield_dGamma += CalculatePlasticStrainRateDerivative(
					MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature)/ CurrentProcessInfo[DELTA_TIME];
				dYield_dGamma += MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] * yield_stress/ MaterialProperties[DENSITY] /
					MaterialProperties[SPECIFIC_HEAT] * CalculateThermalDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
				dYield_dGamma *= GetSqrt23();
				yield_function_gradient = -1.0 * GetSqrt6() * shear_modulus_G - dYield_dGamma;

				// Update gamma and check for convergence
				delta_gamma = -1.0 * yield_function / yield_function_gradient;
				if (std::abs(delta_gamma) < tolerance_delta_gamma) {
					is_converged = true;
					break;
				}
				else gamma += delta_gamma;

				// Bisect increment if out of search bounds
				if (gamma_min > gamma || gamma > gamma_max) gamma = gamma_min + 0.5 * (gamma_max - gamma_min);

				// Update of quantities
				predicted_eps = mEquivalentPlasticStrainOld + GetSqrt23() * gamma; // eps = equivalent plastic strain
				predicted_eps_rate = GetSqrt23() * gamma / CurrentProcessInfo[DELTA_TIME];
				predicted_temperature = mTemperatureOld;
				predicted_temperature += MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] / GetSqrt6() / MaterialProperties[DENSITY]
					/ MaterialProperties[SPECIFIC_HEAT] * (yield_stress + mYieldStressOld) * gamma;

				iteration += 1;
				if (iteration == iteration_limit)
				{
					KRATOS_INFO("Johnson Cook Material Model") << " Johnson Cook iteration limit exceeded\n";
					KRATOS_WATCH(gamma)
					KRATOS_WATCH(delta_gamma)
					KRATOS_WATCH(yield_function)
					KRATOS_ERROR << "Johnson Cook iteration limit exceeded";
				}
			}
			// Correct trial stress
			stress_deviatoric_converged = stress_deviatoric_trial - 2.0 * shear_modulus_G * gamma * flow_direction_normalized;

			// Store plastic deformation quantities
			mEnergyDissipated += gamma / GetSqrt6() / MaterialProperties[DENSITY] * (mYieldStressOld + yield_stress);
			mEquivalentPlasticStrainOld = predicted_eps;
			mPlasticStrainRateOld = predicted_eps_rate;
			mTemperatureOld = predicted_temperature;
			mYieldStressOld = yield_stress;
			mHardeningRatio = yield_stress / mYieldStressVirgin;
		}
		else
		{
			if (yield_stress / mYieldStressVirgin > yield_stress_failure_ratio)
			{
				stress_deviatoric_converged = stress_deviatoric_trial;
			}
			else
			{
				// Particle has failed. It can only take compressive volumetric stresses!
				stress_deviatoric_converged.clear();
				if (stress_hydrostatic_new > 0.0) stress_hydrostatic_new = 0.0;
				mHardeningRatio = 0.0;
			}
		}
		// Update equivalent stress
		Matrix stress_converged = stress_deviatoric_converged + stress_hydrostatic_new * identity;
		mEquivalentStress = std::sqrt(3.0 / 2.0 *
			MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_deviatoric_converged));

		// Store stresses and strains
		MakeStrainStressVectorFromMatrix(stress_converged, StressVector);
		mStrainOld = StrainVector;

		// Udpdate internal energy
		for (size_t i = 0; i < stress_converged.size1(); ++i) {
			for (size_t j = 0; j < stress_converged.size2(); ++j) {
				mEnergyInternal += 0.5 * (stress_converged(i, j) + stress_old(i, j)) * strain_increment(i, j) / MaterialProperties[DENSITY];
			}
		}
	}


	int JohnsonCookThermalPlastic3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const
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


	bool JohnsonCookThermalPlastic3DLaw::CheckParameters(Parameters& rValues)
	{
		return rValues.CheckAllParameters();
	}


	void JohnsonCookThermalPlastic3DLaw::MakeStrainStressVectorFromMatrix(const Matrix& rInput, Vector& rOutput)
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


	void JohnsonCookThermalPlastic3DLaw::MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput)
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


	void JohnsonCookThermalPlastic3DLaw::CheckIsExplicitTimeIntegration(const ProcessInfo& rCurrentProcessInfo)
	{
		const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
			? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
			: false;
		KRATOS_ERROR_IF_NOT(is_explicit) << "The Johnson Cook MPM material law is currently limited to explicit time integration only";
	}


	double JohnsonCookThermalPlastic3DLaw::CalculateHardenedYieldStress(const Properties& MaterialProperties,
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		// Hardening formula is: = (A + B* ep^n) * strain_rate_hardening_factor * thermal_hardening_factor
		double hardened_stress = MaterialProperties[JC_PARAMETER_A] + MaterialProperties[JC_PARAMETER_B] *
			std::pow(EquivalentPlasticStrain, MaterialProperties[JC_PARAMETER_n]);
		hardened_stress *= CalculateStrainRateHardeningFactor(MaterialProperties, PlasticStrainRate);
		hardened_stress *= CalculateThermalHardeningFactor(MaterialProperties, Temperature);

		return hardened_stress;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculateThermalHardeningFactor(const Properties& MaterialProperties, const double Temperature)
	{
		// Calculate thermal hardening factor
		double thermal_hardening_factor;
		if (MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] == 0.0) thermal_hardening_factor = 1.0;
		else if (Temperature < MaterialProperties[REFERENCE_TEMPERATURE]) thermal_hardening_factor = 1.0;
		else if (Temperature >= MaterialProperties[MELD_TEMPERATURE]) thermal_hardening_factor = 0.0;
		else {
			thermal_hardening_factor = 1.0 - std::pow(
			(Temperature - MaterialProperties[REFERENCE_TEMPERATURE]) /
			(MaterialProperties[MELD_TEMPERATURE] - MaterialProperties[REFERENCE_TEMPERATURE])
				, MaterialProperties[JC_PARAMETER_m]);
		}

		return thermal_hardening_factor;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculateStrainRateHardeningFactor(const Properties& MaterialProperties,
		const double PlasticStrainRate)
	{
		// Calculate strain rate hardening factor

		double strain_rate_hardening_factor = 1.0;
		if (PlasticStrainRate > MaterialProperties[REFERENCE_STRAIN_RATE])
		{
			strain_rate_hardening_factor += MaterialProperties[JC_PARAMETER_C] *
				std::log(PlasticStrainRate / MaterialProperties[REFERENCE_STRAIN_RATE]);
		}

		return strain_rate_hardening_factor;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculateThermalDerivative(const Properties& MaterialProperties,
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		double thermal_derivative = 0.0;
		if (MaterialProperties[REFERENCE_TEMPERATURE] <= Temperature
			&& Temperature <= MaterialProperties[MELD_TEMPERATURE]
			&& MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] > 0.0)
		{
			thermal_derivative = -1.0 * MaterialProperties[JC_PARAMETER_m] *
				(MaterialProperties[JC_PARAMETER_A] + MaterialProperties[JC_PARAMETER_B] *
				std::pow(EquivalentPlasticStrain, MaterialProperties[JC_PARAMETER_n])) /
				(Temperature - MaterialProperties[REFERENCE_TEMPERATURE]);
			thermal_derivative *= CalculateStrainRateHardeningFactor(MaterialProperties, PlasticStrainRate);
			thermal_derivative *= std::pow(
			(Temperature - MaterialProperties[REFERENCE_TEMPERATURE]) /
			(MaterialProperties[MELD_TEMPERATURE] - MaterialProperties[REFERENCE_TEMPERATURE]), MaterialProperties[JC_PARAMETER_m]);
		}

		return thermal_derivative;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculatePlasticStrainRateDerivative(const Properties& MaterialProperties,
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		double plastic_strain_rate_derivative = 0.0;

		if (PlasticStrainRate >= MaterialProperties[REFERENCE_STRAIN_RATE])
		{
			plastic_strain_rate_derivative = MaterialProperties[JC_PARAMETER_C] / PlasticStrainRate *
				(MaterialProperties[JC_PARAMETER_A] + MaterialProperties[JC_PARAMETER_B] * std::pow(EquivalentPlasticStrain, MaterialProperties[JC_PARAMETER_n]));
			plastic_strain_rate_derivative *= CalculateThermalHardeningFactor(MaterialProperties, Temperature);
		}

		return plastic_strain_rate_derivative;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculatePlasticStrainDerivative(const Properties& MaterialProperties,
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		double plastic_strain_derivative = MaterialProperties[JC_PARAMETER_n] * MaterialProperties[JC_PARAMETER_B] *
			std::pow(EquivalentPlasticStrain, (MaterialProperties[JC_PARAMETER_n] - 1.0));
		plastic_strain_derivative *= CalculateStrainRateHardeningFactor(MaterialProperties, PlasticStrainRate);
		plastic_strain_derivative *= CalculateThermalHardeningFactor(MaterialProperties, Temperature);

		return plastic_strain_derivative;
	}


	void JohnsonCookThermalPlastic3DLaw::GetLawFeatures(Features& rFeatures)
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


	double& JohnsonCookThermalPlastic3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == MP_TEMPERATURE)
		{
			rValue = mTemperatureOld;
		}
		else if (rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN)
		{
			rValue = mEquivalentPlasticStrainOld;
		}
		else if (rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN_RATE)
		{
			rValue = mPlasticStrainRateOld;
		}
		else if (rThisVariable == MP_HARDENING_RATIO)
		{
			rValue = mHardeningRatio;
		}
		else if (rThisVariable == MP_EQUIVALENT_STRESS)
		{
			rValue = mEquivalentStress;
		}
		else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function GetValue double.";

		return(rValue);
	}


	void JohnsonCookThermalPlastic3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rThisVariable == MP_TEMPERATURE)
		{
			if (mEquivalentPlasticStrainOld > 0)
			{
				KRATOS_ERROR << "Temperature may only be set before any plastic strain occurs.";
			}
			else mTemperatureOld = rValue;
		}
		else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function SetValue double.";
	}
} // Namespace Kratos