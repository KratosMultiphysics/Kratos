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
#include "particle_mechanics_application_variables.h"

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
		if (rThisVariable == TEMPERATURE
			|| rThisVariable == PLASTIC_STRAIN
			|| rThisVariable == PLASTIC_STRAIN_RATE) return true;
		else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function Has double.";

		return false;
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

		mYieldStressOld = CalculateHardenedYieldStress(rMaterialProperties, mEquivalentPlasticStrainOld, mPlasticStrainRateOld, mTemperatureOld);
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

		// Get old stress vector and current strain vector
		const Vector StrainVector = rValues.GetStrainVector();
		Vector StressVector = rValues.GetStressVector();

		// Convert vectors to matrices for easier manipulation
		Matrix stress_old(3, 3);
		Matrix strain_increment(3, 3);
		MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
		MakeStrainStressMatrixFromVector(StressVector, stress_old);

		// Material moduli
		const double shear_modulus_G = MaterialProperties[YOUNG_MODULUS] / (2.0 + 2.0 * MaterialProperties[POISSON_RATIO]);
		const double bulk_modulus_K = MaterialProperties[YOUNG_MODULUS] / (3.0 - 6.0 * MaterialProperties[POISSON_RATIO]);
		const double density = MaterialProperties[DENSITY];

		// Calculate strain increments
		const double strain_increment_trace = strain_increment(0, 0) + strain_increment(1, 1) + strain_increment(2, 2);
		const Matrix strain_increment_hydrostatic = strain_increment_trace / 3.0 * IdentityMatrix(3);
		const Matrix strain_increment_deviatoric = strain_increment - strain_increment_hydrostatic;

		// Calculate current (old) j2 stress
		const double stress_hydrostatic_old = (stress_old(0, 0) + stress_old(1, 1) + stress_old(2, 2)) / 3.0;
		Matrix stress_deviatoric_old = stress_old - stress_hydrostatic_old * IdentityMatrix(3);

		// Calculate trial (predicted) j2 stress
		double stress_hydrostatic_new = stress_hydrostatic_old + bulk_modulus_K * strain_increment_trace / 3.0; // TODO (1) - think this is fixed
		Matrix stress_deviatoric_trial = stress_deviatoric_old + 2.0 * shear_modulus_G * strain_increment_deviatoric;
		const double j2_stress_trial = std::sqrt(3.0 / 2.0 * CalculateMatrixDoubleContraction(stress_deviatoric_trial));

		// Declare deviatoric stress matrix to be used later
		Matrix stress_deviatoric_converged(3, 3);

		// Assume current yield stress is the same as the old (elastic predictor)
		double yield_stress = mYieldStressOld;

		if (j2_stress_trial > yield_stress)
		{
			const double j2_stress_old = std::sqrt(3.0 / 2.0 * CalculateMatrixDoubleContraction(stress_deviatoric_old));

			// Thermal properties
			const double eta = 0.9; // TODO check this
			const double specific_heat_Cp = MaterialProperties[SPECIFIC_HEAT];

			// Newton raphson setup
			double gamma = mGammaOld;
			double gamma_min = 0.0;
			double gamma_max = j2_stress_trial / std::sqrt(6.0) / shear_modulus_G;
			bool is_converged = false;
			const SizeType iteration_limit = 100;
			IndexType iteration = 0;
			const double tolerance = 1e-8;
			double yield_function, yield_function_gradient, dYield_dGamma,
				delta_gamma, predicted_eps, predicted_eps_rate,
				predicted_temperature;

			// Initial prediction of quantities
			predicted_eps = mEquivalentPlasticStrainOld + std::sqrt(2.0 / 3.0) * gamma; // eps = equivalent plastic strain
			predicted_eps_rate = std::sqrt(2.0 / 3.0) * gamma / CurrentProcessInfo[DELTA_TIME];
			predicted_temperature = mTemperatureOld + eta / std::sqrt(6.0) / density / specific_heat_Cp * // TODO check this, [johnson cool umat pdfp169]
				(yield_stress + mYieldStressOld) * gamma;

			// Newton Raphson return mapping loop
			while (!is_converged)
			{
				// Calculate predicted yield stress
				yield_stress = CalculateHardenedYieldStress(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);

				// Compute yield function and derivative
				yield_function = j2_stress_trial - std::sqrt(6.0) * shear_modulus_G * gamma - yield_stress;
				if (yield_function < 0.0) gamma_max = gamma;
				else gamma_min = gamma;
				dYield_dGamma = CalculatePlasticStrainDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
				dYield_dGamma += CalculateThermalDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
				dYield_dGamma += CalculatePlasticStrainRateDerivative(MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature);
				dYield_dGamma *= std::sqrt(2.0 / 3.0);
				yield_function_gradient = -1.0 * std::sqrt(6.0) * shear_modulus_G - dYield_dGamma;

				// Update gamma
				delta_gamma = -1.0 * yield_function / yield_function_gradient;
				gamma += delta_gamma;

				// Bisect increment if out of search bounds
				if (gamma_min > gamma || gamma > gamma_max) gamma = gamma_min + 0.5 * (gamma_max - gamma_min);

				if (std::abs(delta_gamma) < tolerance) {
					is_converged = true;
					break;
				}

				// Update of quantities
				predicted_eps = mEquivalentPlasticStrainOld + std::sqrt(2.0 / 3.0) * gamma; // eps = equivalent plastic strain
				predicted_eps_rate = std::sqrt(2.0 / 3.0) * gamma / CurrentProcessInfo[DELTA_TIME];
				predicted_temperature = mTemperatureOld + eta / std::sqrt(6.0) / density / specific_heat_Cp * // TODO check this, [johnson cool umat pdfp169]
					(yield_stress + mYieldStressOld) * gamma;

				iteration += 1;
				KRATOS_ERROR_IF(iteration == iteration_limit) << "Johnson Cook iteration limit exceeded";
			}

			// Correct trial stress
			Matrix flow_direction_normalized = stress_deviatoric_trial / (j2_stress_trial / std::sqrt(3.0 / 2.0));
			stress_deviatoric_converged = stress_deviatoric_trial - 2.0 * shear_modulus_G * gamma * flow_direction_normalized;

			mEnergyDissipated += gamma / std::sqrt(6.0) / density * (mYieldStressOld + yield_stress);
			mEquivalentPlasticStrainOld = predicted_eps;
			mPlasticStrainRateOld = predicted_eps_rate;
			mTemperatureOld = predicted_temperature;
			mYieldStressOld = yield_stress;
		}
		else
		{
			stress_deviatoric_converged = stress_deviatoric_trial;
		}

		Matrix stress_converged = stress_deviatoric_converged + stress_hydrostatic_new * IdentityMatrix(3);

		// Store stresses and strains
		MakeStrainStressVectorFromMatrix(stress_converged, StressVector);
		mStrainOld = StrainVector;

		// Udpdate internal energy
		for (size_t i = 0; i < stress_converged.size1(); ++i)
		{
			for (size_t j = 0; j < stress_converged.size2(); ++j)
			{
				mEnergyInternal += 0.5 * (stress_converged(i, j) + stress_old(i, j)) * strain_increment(i, j) / density;
			}
		}

		if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
		{
			KRATOS_ERROR << "COMPUTE_CONSTITUTIVE_TENSOR not yet implemented in JohnsonCookThermalPlastic3DLaw";
		}
	}


	int JohnsonCookThermalPlastic3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

		KRATOS_CHECK_VARIABLE_KEY(JC_PARAMETER_A);
		KRATOS_CHECK_VARIABLE_KEY(JC_PARAMETER_B);
		KRATOS_CHECK_VARIABLE_KEY(JC_PARAMETER_m);
		KRATOS_CHECK_VARIABLE_KEY(JC_PARAMETER_n);
		KRATOS_CHECK_VARIABLE_KEY(REFERENCE_STRAIN_RATE);
		KRATOS_CHECK_VARIABLE_KEY(REFERENCE_TEMPERATURE);
		KRATOS_CHECK_VARIABLE_KEY(MELD_TEMPERATURE);

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
		//Constant Parameters of the -- Johnson and Cook --:
		const double A = MaterialProperties[JC_PARAMETER_A];
		const double B = MaterialProperties[JC_PARAMETER_B];
		const double n = MaterialProperties[JC_PARAMETER_n];

		// Hardening formula is: = (A + B* ep^n) * strain_rate_hardening_factor * thermal_hardening_factor
		double hardened_stress = A + B * std::pow(EquivalentPlasticStrain, n);
		hardened_stress *= CalculateStrainRateHardeningFactor(MaterialProperties, PlasticStrainRate);
		hardened_stress *= CalculateThermalHardeningFactor(MaterialProperties, Temperature);

		return hardened_stress;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculateThermalHardeningFactor(const Properties& MaterialProperties, const double Temperature)
	{
		// Calculate thermal hardening factor
		const double m = MaterialProperties[JC_PARAMETER_m];
		const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
		const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];

		double thermal_hardening_factor;
		if (Temperature < ReferenceTemperature) thermal_hardening_factor = 1.0;
		else if (Temperature >= MeldTemperature) thermal_hardening_factor = 0.0;
		else {
			thermal_hardening_factor = 1.0 - std::pow((Temperature - ReferenceTemperature) / (MeldTemperature - ReferenceTemperature), m);
		}

		return thermal_hardening_factor;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculateStrainRateHardeningFactor(const Properties& MaterialProperties, 
		const double PlasticStrainRate)
	{
		// Calculate strain rate hardening factor
		const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];
		const double C = MaterialProperties[JC_PARAMETER_C];

		double strain_rate_hardening_factor = 1.0;
		if (PlasticStrainRate > ReferenceStrainRate)
		{
			strain_rate_hardening_factor += C * std::log(PlasticStrainRate / ReferenceStrainRate);
		}

		return strain_rate_hardening_factor;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculateThermalDerivative(const Properties& MaterialProperties, 
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		const double A = MaterialProperties[JC_PARAMETER_A];
		const double B = MaterialProperties[JC_PARAMETER_B];
		const double C = MaterialProperties[JC_PARAMETER_C];
		const double n = MaterialProperties[JC_PARAMETER_n];
		const double m = MaterialProperties[JC_PARAMETER_m];
		const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
		const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];
		const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];

		double thermal_derivative = 0.0;
		if (ReferenceTemperature <= Temperature && Temperature <= MeldTemperature)
		{
			double strain_rate_hardening_factor = CalculateStrainRateHardeningFactor(MaterialProperties, PlasticStrainRate);
			double thermal_hardening_factor = std::pow((Temperature - ReferenceTemperature) / (MeldTemperature - ReferenceTemperature), m);
			double temp = -1.0 * m * (A + B * std::pow(EquivalentPlasticStrain, n)) / (Temperature - ReferenceTemperature);
			thermal_derivative = temp * strain_rate_hardening_factor * thermal_hardening_factor;
		}

		return thermal_derivative;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculatePlasticStrainRateDerivative(const Properties& MaterialProperties, 
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		const double A = MaterialProperties[JC_PARAMETER_A];
		const double B = MaterialProperties[JC_PARAMETER_B];
		const double C = MaterialProperties[JC_PARAMETER_C];
		const double n = MaterialProperties[JC_PARAMETER_n];
		const double m = MaterialProperties[JC_PARAMETER_m];
		const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
		const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];
		const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];

		double plastic_strain_rate_derivative = 0.0;

		if (PlasticStrainRate >= ReferenceStrainRate)
		{
			double thermal_hardening_factor = CalculateThermalHardeningFactor(MaterialProperties, Temperature);

			plastic_strain_rate_derivative = C / PlasticStrainRate * (A + B * std::pow(EquivalentPlasticStrain, n)) * thermal_hardening_factor;
		}

		return plastic_strain_rate_derivative;
	}


	double JohnsonCookThermalPlastic3DLaw::CalculatePlasticStrainDerivative(const Properties& MaterialProperties, 
		const double EquivalentPlasticStrain, const double PlasticStrainRate, const double Temperature)
	{
		const double A = MaterialProperties[JC_PARAMETER_A];
		const double B = MaterialProperties[JC_PARAMETER_B];
		const double C = MaterialProperties[JC_PARAMETER_C];
		const double n = MaterialProperties[JC_PARAMETER_n];
		const double m = MaterialProperties[JC_PARAMETER_m];
		const double ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
		const double MeldTemperature = MaterialProperties[MELD_TEMPERATURE];
		const double ReferenceStrainRate = MaterialProperties[REFERENCE_STRAIN_RATE];

		double plastic_strain_derivative = 0.0;

		// Calculate thermal hardening factor
		double thermal_hardening_factor = CalculateThermalHardeningFactor(MaterialProperties, Temperature);

		// Calculate strain rate hardening factor
		double strain_rate_hardening_factor = 1.0;
		if (PlasticStrainRate > ReferenceStrainRate)
		{
			strain_rate_hardening_factor += C * std::log(PlasticStrainRate / ReferenceStrainRate);
		}

		plastic_strain_derivative = n * B * std::pow(EquivalentPlasticStrain, (n - 1.0)) * strain_rate_hardening_factor * thermal_hardening_factor;

		return plastic_strain_derivative;
	}


	void JohnsonCookThermalPlastic3DLaw::GetLawFeatures(Features& rFeatures)
	{
		//Set the type of law
		rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
		rFeatures.mOptions.Set(FINITE_STRAINS);
		rFeatures.mOptions.Set(ISOTROPIC);

		//Set strain measure required by the consitutive law
		rFeatures.mStrainMeasures.push_back(StrainMeasure_GreenLagrange);

		//Set the strain size
		rFeatures.mStrainSize = GetStrainSize();

		//Set the spacedimension
		rFeatures.mSpaceDimension = WorkingSpaceDimension();
	}


	double& JohnsonCookThermalPlastic3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == TEMPERATURE)
		{
			rValue = mTemperatureOld;
		}
		else if (rThisVariable == PLASTIC_STRAIN)
		{
			rValue = mEquivalentPlasticStrainOld;
		}
		else if (rThisVariable == PLASTIC_STRAIN_RATE)
		{
			rValue = mEquivalentPlasticStrainOld;
		}
		KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in Johnson Cook 3D material law function GetValue double.";

		return(rValue);
	}


	void JohnsonCookThermalPlastic3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rThisVariable == TEMPERATURE)
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