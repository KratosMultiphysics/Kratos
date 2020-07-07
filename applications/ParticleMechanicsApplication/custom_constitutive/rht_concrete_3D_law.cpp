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
			|| rThisVariable == MP_EQUIVALENT_STRESS
			|| rThisVariable == MP_DAMAGE
			|| rThisVariable == MP_COMPACTION_RATIO)
			return true;
		else return false;
	}


	void RHTConcrete3DLaw::InitializeMaterial(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
	{
		BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

		mStrainOld = ZeroVector(GetStrainSize());
		mDamage = 0.0;
		mAlpha = 1.1;
		mEquivalentStress = 0.0;
		mEquivalentPlasticStrain = 0.0;
		mEquivalentPlasticStrainRate = 0.0;
		mDensityInitial = rMaterialProperties[DENSITY];
	}


	void RHTConcrete3DLaw::CalculateMaterialResponseKirchhoff(Kratos::ConstitutiveLaw::Parameters& rValues)
	{
		// Check if the constitutive parameters are passed correctly to the law calculation
		CheckParameters(rValues);

		// Get Values to compute the constitutive law:
		Flags& Options = rValues.GetOptions();
		KRATOS_ERROR_IF(Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
			<< "The RHTConcrete3DLaw cannot accept the deformation gradient F as a strain input."
			<< " Please set the strain vector and set USE_ELEMENT_PROVIDED_STRAIN = False in the element.";
		const ProcessInfo& process_info = rValues.GetProcessInfo();
		CheckIsExplicitTimeIntegration(process_info);
		const Properties& mat_props = rValues.GetMaterialProperties();

		// Get old stress vector and current strain vector
		const Vector StrainVector = rValues.GetStrainVector();
		Vector& StressVector = rValues.GetStressVector();
		const Matrix identity = (GetStrainSize() == 3) ? IdentityMatrix(2) : IdentityMatrix(3);

		// Convert vectors to matrices for easier manipulation
		Matrix stress = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		Matrix strain_increment = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		Matrix strain_new = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		MPMStressPrincipalInvariantsUtility::MakeStrainStressMatrixFromVector(
			(StrainVector - mStrainOld), strain_increment,GetStrainSize());
		MPMStressPrincipalInvariantsUtility::MakeStrainStressMatrixFromVector(
			StressVector, stress, GetStrainSize());
		MPMStressPrincipalInvariantsUtility::MakeStrainStressMatrixFromVector(
			StrainVector, strain_new, GetStrainSize());

		// Calculate deviatoric quantities
		const Matrix strain_increment_deviatoric = strain_increment -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(strain_increment) / 3.0 * identity;
		const Matrix stress_deviatoric_old = stress -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress) / 3.0 * identity;
		Matrix stress_deviatoric = stress_deviatoric_old + 2.0 * mat_props[SHEAR_MODULUS] * strain_increment_deviatoric;

		// Calculate pressure from Equation of State
		const double pressure = CalculatePressureFromEOS(mat_props, strain_new, stress);

		// Prepare elastic predictor
		const double eff_stress_trial = std::sqrt(3.0 / 2.0 *
			MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_deviatoric));
		const double lode_angle =
			MPMStressPrincipalInvariantsUtility::CalculateLodeAngleFromDeviatoricStressTensor(stress_deviatoric);

		// Assume elastic, EPSrate = 0.0
		double eps_rate = 0.0;
		double limit_elastic = CalculateElasticLimit(pressure, lode_angle,
			CalculateRateFactors(eps_rate, mat_props), mEquivalentPlasticStrain, mat_props);

		if (eff_stress_trial > limit_elastic)
		{
			// Plastic deformation, need to correct

			double delta_eps_min = 0.0;
			double delta_eps_max = eff_stress_trial/3.0* mat_props[SHEAR_MODULUS];
			double delta_eps = 0.5 * (delta_eps_min + delta_eps_max);
			eps_rate = delta_eps / process_info[DELTA_TIME];
			double eps_trial = mEquivalentPlasticStrain + delta_eps;

			const SizeType iteration_limit = 100;
			IndexType iteration = 1;
			bool is_converged = false;

			const Matrix flow_direction = stress_deviatoric /
				(eff_stress_trial / std::sqrt(3.0 / 2.0));

			double limit_failure, shear_mod_plastic, eps_harden_limit, eps_hard_ratio,
				limit_current_surface;
			double delta_damage = 0.0;
			double damage_trial = mDamage;

			while (!is_converged)
			{
				const array_1d<double, 2> rate_factors = CalculateRateFactors(eps_rate, mat_props);

				limit_elastic = CalculateElasticLimit(pressure, lode_angle, rate_factors, eps_trial, mat_props);
				limit_failure = CalculateFailureLimit(pressure, lode_angle, rate_factors, mat_props);

				shear_mod_plastic = mat_props[SHEAR_MODULUS] * mat_props[RHT_SHEAR_MOD_REDUCTION_FACTOR];

				// Equivalent plastic strain at fully hardened state ref[3] eqn5.32
				eps_harden_limit = (limit_failure - limit_elastic) / 3.0 / shear_mod_plastic;
				eps_hard_ratio = std::min(eps_trial / eps_harden_limit, 1.0);

				// Check if we have reached the hardening limit, damage accumulates after this
				if (eps_hard_ratio == 1.0)
				{
					const double eps_fail = CalculateFailureStrain(pressure,
						damage_trial, rate_factors, mat_props);
					// Calculate damage increment, ref[1]eqn21, ref[2]
					damage_trial = (eps_trial - eps_harden_limit) / eps_fail;
					if (damage_trial > 1.0) damage_trial = 1.0;
					if (damage_trial < mDamage) damage_trial = mDamage; // Reversal of damage not allowed!

					const double limit_residual = CalculateResidualLimit(pressure, mat_props);
					// Our current surface lies in between the failure and residual surface.
					// Implemented as per ref[4] and ref[2]
					limit_current_surface = damage_trial * limit_residual + (1.0 - damage_trial) * limit_failure;
				}
				else
				{
					// Our current surface lies in between the elastic and failure surface.
					// Implemented as per ref[4] and ref[2]
					limit_current_surface = eps_hard_ratio * limit_failure + (1.0 - eps_hard_ratio) * limit_elastic;
				}

				const double consistency_condition = eff_stress_trial -
					limit_current_surface - 3.0 * mat_props[SHEAR_MODULUS] * delta_eps;
				if (consistency_condition < 0.0) delta_eps_max = delta_eps;
				else delta_eps_min = delta_eps;
				const double delta_eps_increment = 0.5 * (delta_eps_max - delta_eps_min);

				if (std::abs(delta_eps_increment) < 1e-12)
				{
					if (std::abs(consistency_condition)< mat_props[RHT_COMPRESSIVE_STRENGTH]/1e6)
					{
						is_converged = true;

						break;
					}
				}

				delta_eps = delta_eps_min + delta_eps_increment;
				eps_trial = mEquivalentPlasticStrain + delta_eps;
				eps_rate = delta_eps / process_info[DELTA_TIME];

				iteration += 1;
				if (iteration > iteration_limit)
				{
					KRATOS_INFO("RHTConcrete3DLaw::CalculateMaterialResponseKirchhoff")
						<< "Iteration = " << iteration << "\n"
						<< "delta_eps_increment = " << delta_eps_increment << "\n"
						<< "damage_trial = " << damage_trial << "\n";
					KRATOS_ERROR << "Iteration limit exceeded\n";
				}
			}

			// Update strength quantities
			mEquivalentPlasticStrain = eps_trial;
			mDamage = damage_trial;
			stress_deviatoric -= std::sqrt(6.0) * mat_props[SHEAR_MODULUS] * delta_eps * flow_direction;
		}

		// Store variables
		mEquivalentStress = std::sqrt(3.0 / 2.0 *
			MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_deviatoric));
		mEquivalentPlasticStrainRate = eps_rate;

		// Recombine stress matrix
		stress = stress_deviatoric - pressure * identity;

		// Store stresses and strains
		MPMStressPrincipalInvariantsUtility::MakeStrainStressVectorFromMatrix(
			stress, StressVector,GetStrainSize());
		mStrainOld = StrainVector;
	}


	int RHTConcrete3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

		// Note, expected values taken from ref[1] for standard (C30/37) concrete with compressive strength of 35MPa

		KRATOS_ERROR_IF (SHEAR_MODULUS.Key()==0 || rMaterialProperties[SHEAR_MODULUS] < 0.0)
			<< "SHEAR_MODULUS has key zero or invalid value (expected positive number ~16700MPa)" << std::endl;

		KRATOS_ERROR_IF (DENSITY.Key()==0 || rMaterialProperties[DENSITY] < 0.0)
			<< "DENSITY has key zero or invalid value (expected positive number ~2314kg/m3)" << std::endl;

		KRATOS_ERROR_IF (RHT_A.Key()==0 || rMaterialProperties[RHT_A] < 0.0)
			<< "RHT_A has key zero or invalid value (expected positive number ~1.6)" << std::endl;

		KRATOS_ERROR_IF (RHT_N.Key()==0 || rMaterialProperties[RHT_N] < 0.0)
			<< "RHT_N has key zero or invalid value (expected positive number ~0.61)" << std::endl;

		KRATOS_ERROR_IF (RHT_COMPRESSIVE_STRENGTH.Key()==0 || rMaterialProperties[RHT_COMPRESSIVE_STRENGTH] < 0.0)
			<< "RHT_COMPRESSIVE_STRENGTH has key zero or invalid value (expected positive number ~35MPa)" << std::endl;

		KRATOS_ERROR_IF (RHT_RELATIVE_SHEAR_STRENGTH.Key()==0 || rMaterialProperties[RHT_RELATIVE_SHEAR_STRENGTH] < 0.0)
			<< "RHT_RELATIVE_SHEAR_STRENGTH has key zero or invalid value (expected positive number ~0.18)" << std::endl;

		KRATOS_ERROR_IF (RHT_RELATIVE_TENSILE_STRENGTH.Key()==0 || rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH] < 0.0)
			<< "RHT_RELATIVE_TENSILE_STRENGTH has key zero or invalid value (expected positive number ~0.10)" << std::endl;

		KRATOS_ERROR_IF (RHT_Q0.Key()==0 || rMaterialProperties[RHT_Q0] < 0.0)
			<< "RHT_Q0 has key zero or invalid value (expected positive number ~0.6805)" << std::endl;

		KRATOS_ERROR_IF (RHT_B.Key()==0 || rMaterialProperties[RHT_B] < 0.0)
			<< "RHT_B has key zero or invalid value (expected positive number ~0.0105)" << std::endl;

		KRATOS_ERROR_IF (REFERENCE_TENSION_STRAIN_RATE.Key()==0 || rMaterialProperties[REFERENCE_TENSION_STRAIN_RATE] < 0.0)
			<< "REFERENCE_TENSION_STRAIN_RATE has key zero or invalid value (expected positive number ~3e-6)" << std::endl;

		KRATOS_ERROR_IF (REFERENCE_COMPRESSION_STRAIN_RATE.Key()==0 || rMaterialProperties[REFERENCE_COMPRESSION_STRAIN_RATE] < 0.0)
			<< "REFERENCE_COMPRESSION_STRAIN_RATE has key zero or invalid value (expected positive number ~3e-5)" << std::endl;

		KRATOS_ERROR_IF (RHT_GC_STAR.Key()==0 || rMaterialProperties[RHT_GC_STAR] < 0.0)
			<< "RHT_GC_STAR has key zero or invalid value (expected positive number ~0.53)" << std::endl;

		KRATOS_ERROR_IF (RHT_GT_STAR.Key()==0 || rMaterialProperties[RHT_GT_STAR] < 0.0)
			<< "RHT_GT_STAR has key zero or invalid value (expected positive number ~0.7)" << std::endl;

		KRATOS_ERROR_IF (RHT_SHEAR_MOD_REDUCTION_FACTOR.Key()==0 || rMaterialProperties[RHT_SHEAR_MOD_REDUCTION_FACTOR] < 0.0)
			<< "RHT_SHEAR_MOD_REDUCTION_FACTOR has key zero or invalid value (expected positive number ~0.5)" << std::endl;

		KRATOS_ERROR_IF (RHT_D1.Key()==0 || rMaterialProperties[RHT_D1] < 0.0)
			<< "RHT_D1 has key zero or invalid value (expected positive number ~0.04)" << std::endl;

		KRATOS_ERROR_IF (RHT_D2.Key()==0 || rMaterialProperties[RHT_D2] < 0.0)
			<< "RHT_D2 has key zero or invalid value (expected positive number ~1.0)" << std::endl;

		KRATOS_ERROR_IF (RHT_MIN_DAMAGED_RESIDUAL_STRAIN.Key()==0 || rMaterialProperties[RHT_MIN_DAMAGED_RESIDUAL_STRAIN] < 0.0)
			<< "RHT_MIN_DAMAGED_RESIDUAL_STRAIN has key zero or invalid value (expected positive number ~0.01)" << std::endl;

		KRATOS_ERROR_IF (RHT_AF.Key()==0 || rMaterialProperties[RHT_AF] < 0.0)
			<< "RHT_AF has key zero or invalid value (expected positive number ~1.6)" << std::endl;

		KRATOS_ERROR_IF (RHT_NF.Key()==0 || rMaterialProperties[RHT_NF] < 0.0)
			<< "RHT_NF has key zero or invalid value (expected positive number ~0.61)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_A1.Key()==0 || rMaterialProperties[RHT_EOS_A1] < 0.0)
			<< "RHT_EOS_A1 has key zero or invalid value (expected positive number ~3.527e10)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_A2.Key()==0 || rMaterialProperties[RHT_EOS_A2] < 0.0)
			<< "RHT_EOS_A2 has key zero or invalid value (expected positive number ~3.958e10)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_A3.Key()==0 || rMaterialProperties[RHT_EOS_A3] < 0.0)
			<< "RHT_EOS_A3 has key zero or invalid value (expected positive number ~9.04e9)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_B0.Key()==0 || rMaterialProperties[RHT_EOS_B0] < 0.0)
			<< "RHT_EOS_B0 has key zero or invalid value (expected positive number ~1.22)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_B1.Key()==0 || rMaterialProperties[RHT_EOS_B1] < 0.0)
			<< "RHT_EOS_B1 has key zero or invalid value (expected positive number ~1.22)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_T1.Key()==0 || rMaterialProperties[RHT_EOS_T1] < 0.0)
			<< "RHT_EOS_T1 has key zero or invalid value (expected positive number ~3.527e10)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_T2.Key()==0 || rMaterialProperties[RHT_EOS_T2] < 0.0)
			<< "RHT_EOS_T2 has key zero or invalid value (expected positive or zero number ~0.0)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_ALPHA0.Key()==0 || rMaterialProperties[RHT_EOS_ALPHA0] < 1.0)
			<< "RHT_EOS_ALPHA0 has key zero or invalid value (expected positive number > 1.0 ~1.1884)" << std::endl;

		KRATOS_ERROR_IF (RHT_EOS_NP.Key()==0 || rMaterialProperties[RHT_EOS_NP] < 0.0)
			<< "RHT_EOS_NP has key zero or invalid value (expected positive number ~3.0)" << std::endl;

		KRATOS_ERROR_IF (RHT_CRUSH_PRESSURE.Key()==0 || rMaterialProperties[RHT_CRUSH_PRESSURE] < 0.0)
			<< "RHT_CRUSH_PRESSURE has key zero or invalid value (expected positive number ~33MPa)" << std::endl;

		KRATOS_ERROR_IF (RHT_COMPACTION_PRESSURE.Key()==0 || rMaterialProperties[RHT_COMPACTION_PRESSURE] < 0.0)
			<< "RHT_COMPACTION_PRESSURE has key zero or invalid value (expected positive number ~6000MPa)" << std::endl;

		if (check_base > 1) return 1;
		return 0;
	}


	bool RHTConcrete3DLaw::CheckParameters(Parameters& rValues)
	{
		return rValues.CheckAllParameters();
	}


	void RHTConcrete3DLaw::CheckIsExplicitTimeIntegration(const ProcessInfo& rCurrentProcessInfo)
	{
		const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
			? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
			: false;
		KRATOS_ERROR_IF_NOT(is_explicit) << "The Johnson Cook MPM material law is currently limited to explicit time integration only";
	}


	double RHTConcrete3DLaw::CalculatePressureFromEOS(const Properties& rMaterialProperties,
		const Matrix& rStrain, const Matrix& rStress)
	{
		// Compute pressure from RHT p-alpha EOS. Ref [2] (modified)

		const double specific_energy_deviatoric =
			CalculateDeviatoricSpecificEnergy(rMaterialProperties, rStrain, rStress);
		const double strain_hydrostatic = MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(rStrain) / 3.0;

		const SizeType iteration_limit = 100;
		IndexType iteration = 1;
		bool is_converged = false;
		bool alpha_smoothing = true; // used to prevent stick-slipping in solving alpha
		double pressure = -1.0 * MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(rStress) / 3.0;
		double pressure_old = pressure;
		double alpha_trial = 0.0;
		double alpha_old = (rMaterialProperties[RHT_EOS_ALPHA0] - 1.0) / 2.0;

		while (!is_converged)
		{
			// Compute compaction parameter alpha
			alpha_trial = 1.0 + (rMaterialProperties[RHT_EOS_ALPHA0] - 1.0)*
				std::pow((
				(rMaterialProperties[RHT_COMPACTION_PRESSURE]-pressure)/
				(rMaterialProperties[RHT_COMPACTION_PRESSURE]-rMaterialProperties[RHT_CRUSH_PRESSURE])),
				rMaterialProperties[RHT_EOS_NP]);
			alpha_trial = std::min(alpha_trial, mAlpha);
			alpha_trial = std::min(alpha_trial, rMaterialProperties[RHT_EOS_ALPHA0]);
			alpha_trial = std::max(1.0, alpha_trial);
			if (alpha_smoothing) alpha_trial = 0.5* (alpha_trial + alpha_old);

			const double nu = alpha_trial * rMaterialProperties[DENSITY] /
				rMaterialProperties[RHT_EOS_ALPHA0] / mDensityInitial - 1.0;

			if (nu >= 0.0) {
				// Compression
				pressure = (rMaterialProperties[RHT_EOS_B0] + nu * rMaterialProperties[RHT_EOS_B1]) * alpha_trial *
					rMaterialProperties[DENSITY] * specific_energy_deviatoric + nu * rMaterialProperties[RHT_EOS_A1]
					+ nu * nu * rMaterialProperties[RHT_EOS_A2] + nu * nu * nu * rMaterialProperties[RHT_EOS_A3];
				pressure /= (1.0 - 1.5 * strain_hydrostatic * alpha_trial *
					(rMaterialProperties[RHT_EOS_B0] + nu * rMaterialProperties[RHT_EOS_B1]));
			}
			else {
				// Tension
				pressure = rMaterialProperties[RHT_EOS_B0] * alpha_trial * rMaterialProperties[DENSITY] *
					specific_energy_deviatoric + nu * rMaterialProperties[RHT_EOS_T1]
					+ nu * nu * rMaterialProperties[RHT_EOS_T2];
				pressure /= (1.0 - 1.5 * strain_hydrostatic * alpha_trial *rMaterialProperties[RHT_EOS_B0]);
			}
			pressure /= alpha_trial;

			if (std::abs(pressure-pressure_old) < 1e-3)
			{
				is_converged = true;
				break;
			}

			iteration += 1;
			if (iteration > iteration_limit)
			{
				KRATOS_INFO("CalculatePressureFromEOS") << "Iteration = " << iteration << "\n"
					<< "Pressure = " << pressure << "\n"
					<< "Pressure_old = " << pressure_old << "\n"
					<< "Alpha_trial = " << alpha_trial << "\n"
					<< "Alpha_old = " << alpha_old << "\n"
					<< "Error = " << (pressure - pressure_old) << "\n";
				KRATOS_ERROR << "CalculatePressureFromEOS iteration limit exceeded\n";
			}
			pressure_old = pressure;
			alpha_old = alpha_trial;
		}

		mAlpha = alpha_trial;

		// Update cap pressure, ref [2]
		mPoreCrushPressure = rMaterialProperties[RHT_COMPACTION_PRESSURE] -
			(rMaterialProperties[RHT_COMPACTION_PRESSURE] - rMaterialProperties[RHT_CRUSH_PRESSURE]) *
			std::pow(
			(mAlpha - 1.0) / (rMaterialProperties[RHT_EOS_ALPHA0] - 1.0),
				(1.0 / rMaterialProperties[RHT_EOS_NP]));

		return pressure;
	}

	double RHTConcrete3DLaw::CalculateDeviatoricSpecificEnergy(const Properties& rMaterialProperties,
		const Matrix& rStrain, const Matrix& rStress)
	{
		// Specific internal energy may be split into deviatoric and hydrostatic components
		// https://www.continuummechanics.org/vonmisesstress.html

		const Matrix stress_dev = rStress -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(rStress) / 3.0 * IdentityMatrix(rStress.size1());
		const Matrix strain_dev = rStrain -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(rStrain) / 3.0 * IdentityMatrix(rStrain.size1());

		double deviatoric_specific_energy = 0.0;

		for (size_t i = 0; i < rStress.size1(); ++i)
		{
			for (size_t j = 0; j < rStress.size2(); ++j)
			{
				deviatoric_specific_energy += stress_dev(i, j) * strain_dev(i, j);
			}
		}
		deviatoric_specific_energy /= rMaterialProperties[DENSITY];

		return deviatoric_specific_energy;
	}

	double RHTConcrete3DLaw::CalculateElasticLimit(const double Pressure, const double LodeAngle,
		const array_1d<double, 2>& RateFactors, const double EPS, const Properties& rMatProps)
	{
		// ref[1] eqn 14
		const double pressure_star = Pressure / rMatProps[RHT_COMPRESSIVE_STRENGTH];

		// Calculate factors
		const double fe_elastic_factor = CalculateFeElasticFactor(pressure_star,
			RateFactors, rMatProps);

		const double fc_cap_factor = CalculateFcCapFactor(pressure_star,
			EPS, RateFactors, rMatProps);

		// Rate factor uses p*, not p*/Fe as per ref[1] eqn14
		const double fr_rate_factor = CalculateFrRateFactor(pressure_star,
			RateFactors, rMatProps);

		const double normalized_yield = CalculateNormalizedYield(pressure_star/ fe_elastic_factor,
			fr_rate_factor, rMatProps);

		const double q_deviatoric_shape_factor = CalculateDeviatoricShapeFactorQ(pressure_star, rMatProps);
		const double r_triaxiality = CalculateTriaxialityR(LodeAngle, q_deviatoric_shape_factor);

		const double elastic_limit = rMatProps[RHT_COMPRESSIVE_STRENGTH] * normalized_yield *
			r_triaxiality * fe_elastic_factor * fc_cap_factor;

		return elastic_limit;
	}

	double RHTConcrete3DLaw::CalculateFailureLimit(const double Pressure, const double LodeAngle,
		const array_1d<double, 2>& RateFactors, const Properties& rMaterialProperties)
	{
		// ref[1] eqn5
		const double pressure_star = Pressure / rMaterialProperties[RHT_COMPRESSIVE_STRENGTH];

		const double fr_rate_factor = CalculateFrRateFactor(pressure_star, RateFactors,
			rMaterialProperties);

		const double normalized_yield = CalculateNormalizedYield(pressure_star,
			fr_rate_factor, rMaterialProperties);

		const double q_deviatoric_shape_factor = CalculateDeviatoricShapeFactorQ(pressure_star, rMaterialProperties);
		const double r_triaxiality = CalculateTriaxialityR(LodeAngle, q_deviatoric_shape_factor);

		const double failure_limit = rMaterialProperties[RHT_COMPRESSIVE_STRENGTH] * normalized_yield *
			r_triaxiality;

		return failure_limit;
	}

	double RHTConcrete3DLaw::CalculateResidualLimit(const double Pressure, const Properties& rMatProps)
	{
		//ref[1] eqn20
		if (Pressure <= 0.0) return 0.0;
		else return rMatProps[RHT_AF] *
				std::pow(Pressure / rMatProps[RHT_COMPRESSIVE_STRENGTH], rMatProps[RHT_NF]);
	}

	const double RHTConcrete3DLaw::CalculateFeElasticFactor(const double PressureStar,
		const array_1d<double, 2> RateFactors, const Properties& rMaterialProperties)
	{
		// ref[1] eqn15
		const double gc_star = rMaterialProperties[RHT_GC_STAR];
		const double gt_star = rMaterialProperties[RHT_GT_STAR];
		const double ft_star = rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH];

		double fe_elastic_factor = 0.0;

		if (3.0* PressureStar >= RateFactors[1] * gc_star)
		{
			fe_elastic_factor = gc_star;
		}
		else if (3.0 * PressureStar >= -1.0* RateFactors[0] * gt_star * ft_star)
		{
			fe_elastic_factor = gc_star - (gt_star - gc_star) * (3.0 * PressureStar - RateFactors[1] * gc_star)/
				(RateFactors[1] * gc_star + RateFactors[0]  * gt_star * ft_star);
		}
		else
		{
			fe_elastic_factor = gt_star;
		}

		return fe_elastic_factor;
	}

	const double RHTConcrete3DLaw::CalculateFcCapFactor(const double PressureStar,
		const double EPS, const array_1d<double, 2> RateFactors,
		const Properties& rMaterialProperties)
	{
		//ref [1] eqn16
		const double degraded_shear_mod = rMaterialProperties[SHEAR_MODULUS] *
			rMaterialProperties[RHT_SHEAR_MOD_REDUCTION_FACTOR];
		const double initial_pressure_star =
			RateFactors[1] * rMaterialProperties[RHT_GC_STAR] / 3.0 +
			degraded_shear_mod * EPS / rMaterialProperties[RHT_COMPRESSIVE_STRENGTH];
		const double pore_crush_star = mPoreCrushPressure / rMaterialProperties[RHT_COMPRESSIVE_STRENGTH];

		double fc_cap_factor;
		if (PressureStar >= pore_crush_star || std::abs(PressureStar - pore_crush_star) < 1e-6)
			fc_cap_factor = 0.0;
		else if (PressureStar >= initial_pressure_star) {
			fc_cap_factor = std::sqrt(1.0 - std::pow((PressureStar - initial_pressure_star) /
			(pore_crush_star - initial_pressure_star), 2.0));
		}
		else fc_cap_factor = 1.0;

		return fc_cap_factor;
	}

	const double RHTConcrete3DLaw::CalculateFrRateFactor(const double PressureStarForRate,
		const array_1d<double, 2> RateFactors, const Properties& rMaterialProperties)
	{
		// ref[1] eqn10

		double fr_rate_factor;
		if (3.0 * PressureStarForRate >= RateFactors[1]) fr_rate_factor = RateFactors[1];
		else if (3.0 * PressureStarForRate >= -1.0 * RateFactors[0] * rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH])
		{
			fr_rate_factor = RateFactors[1] - (RateFactors[0] - RateFactors[1]) *
				(3.0 * PressureStarForRate - RateFactors[1]) /
				(RateFactors[1] + RateFactors[0] * rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH]);
		}
		else fr_rate_factor = RateFactors[0];

		return fr_rate_factor;
	}

	const double RHTConcrete3DLaw::CalculateNormalizedYield(const double PressureStar,
		const double RateFactor, const Properties& rMaterialProperties)
	{
		const double A = rMaterialProperties[RHT_A];
		const double n = rMaterialProperties[RHT_N];
		const double ft_star = rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH];
		const double fs_star = rMaterialProperties[RHT_RELATIVE_SHEAR_STRENGTH];

		const array_1d<double, 2> q_triaxiality_factors =
			CalculateTriaxialityQs(PressureStar, rMaterialProperties);
		const double HTL =
			GetHugoniotTensileLimit(RateFactor, q_triaxiality_factors, rMaterialProperties);

		double normalized_yield;
		if (3.0 * PressureStar >= RateFactor)
		{
			normalized_yield = A * std::pow(
			(PressureStar - RateFactor/3.0+
				std::pow((A/RateFactor),(-1.0/n))
				), n);
		}
		else if (3.0 * PressureStar >= 0.0)
		{
			normalized_yield = RateFactor * fs_star / q_triaxiality_factors[0]
				+ 3.0 * PressureStar * (1.0 - fs_star / q_triaxiality_factors[0]);
		}
		else if (PressureStar >= HTL)
		{
			normalized_yield = RateFactor * fs_star / q_triaxiality_factors[0]
				- 3.0 * PressureStar *
				(1.0/ q_triaxiality_factors[1] - fs_star / q_triaxiality_factors[0]/ft_star);
		}
		else normalized_yield = 0.0;

		return normalized_yield;
	}

	const array_1d<double, 2> RHTConcrete3DLaw::CalculateRateFactors(const double EPSrate,
		const Properties& rMaterialProperties)
	{
		// ref[1] eqn11
		array_1d<double, 2> rate_factors;
		rate_factors[0] = 1.0; // tension
		rate_factors[1] = 1.0; // compression
		if (EPSrate > 0.0) {
			const double beta_t = 1.0 / (10.0 + 0.5 * rMaterialProperties[RHT_COMPRESSIVE_STRENGTH] / 1e6);
			const double beta_c = 1.0 / (5.0 + 0.75 * rMaterialProperties[RHT_COMPRESSIVE_STRENGTH] / 1e6);

			rate_factors[0] = std::pow(EPSrate/ rMaterialProperties[REFERENCE_TENSION_STRAIN_RATE], beta_t);
			rate_factors[1] = std::pow(EPSrate/ rMaterialProperties[REFERENCE_COMPRESSION_STRAIN_RATE], beta_c);
		}
		return rate_factors;
	}

	const array_1d<double, 2> RHTConcrete3DLaw::CalculateTriaxialityQs(const double PressureStar,
		const Properties& rMaterialProperties)
	{
		array_1d<double, 2> q_triaxiality_factors;
		const double q_deviatoric_shape_factor =
			CalculateDeviatoricShapeFactorQ(PressureStar, rMaterialProperties);

		q_triaxiality_factors[0] = CalculateTriaxialityR(Globals::Pi / 6.0, q_deviatoric_shape_factor);
		q_triaxiality_factors[1] = q_deviatoric_shape_factor;

		return q_triaxiality_factors;
	}

	const double RHTConcrete3DLaw::CalculateDeviatoricShapeFactorQ(const double PressureStar,
		const Properties& rMaterialProperties)
	{
		double q_deviatoric_shape_factor = rMaterialProperties[RHT_Q0] +
			rMaterialProperties[RHT_B] * PressureStar;
		if (q_deviatoric_shape_factor < 0.5) q_deviatoric_shape_factor = 0.5;
		if (q_deviatoric_shape_factor > 1.0) q_deviatoric_shape_factor = 1.0;

		return q_deviatoric_shape_factor;
	}

	const double RHTConcrete3DLaw::CalculateTriaxialityR(const double LodeAngle, const double TriaxialityQ)
	{
		// ref[1] eqn8
		double r_triaxiality;

		r_triaxiality = 2.0 * (1.0 - TriaxialityQ * TriaxialityQ) * std::cos(LodeAngle) + (2.0 * TriaxialityQ - 1.0) *
			std::sqrt(4.0 * (1.0 - TriaxialityQ * TriaxialityQ) * std::cos(LodeAngle) * std::cos(LodeAngle)
				+ 5.0 * TriaxialityQ * TriaxialityQ - 4.0 * TriaxialityQ);

		r_triaxiality /= (4.0 * (1.0 - TriaxialityQ * TriaxialityQ) * std::cos(LodeAngle) * std::cos(LodeAngle) +
			(1.0 - 2.0 * TriaxialityQ) * (1.0 - 2.0 * TriaxialityQ));

		return r_triaxiality;
	}

	const double RHTConcrete3DLaw::GetHugoniotTensileLimit(const double RateFactor,
		const array_1d<double, 2> TriaxialityQs, const Properties& rMaterialProperties)
	{
		// ref [1] eqn7
		const double HTL = RateFactor * TriaxialityQs[1] * rMaterialProperties[RHT_RELATIVE_SHEAR_STRENGTH] *
			rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH] / 3.0 /
			(TriaxialityQs[0]* rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH] -
				TriaxialityQs[1] * rMaterialProperties[RHT_RELATIVE_SHEAR_STRENGTH]);

		return HTL;
	}

	const double RHTConcrete3DLaw::CalculateFailureStrain(const double Pressure,
		const double Damage, const array_1d<double, 2>& RateFactors, const Properties& rMatProps)
	{
		// ref [1] eqn22
		const double pressure_star = Pressure / rMatProps[RHT_COMPRESSIVE_STRENGTH];
		const double fr_rate_factor = CalculateFrRateFactor(pressure_star, RateFactors, rMatProps);
		const array_1d<double, 2> q_triaxiality_factors =
			CalculateTriaxialityQs(pressure_star, rMatProps);
		const double HTL =
			GetHugoniotTensileLimit(pressure_star, q_triaxiality_factors, rMatProps);

		double eps_failure = rMatProps[RHT_D1] *
			std::pow((pressure_star - (1.0 - Damage) * HTL), rMatProps[RHT_D2]);
		eps_failure = std::max(eps_failure, rMatProps[RHT_MIN_DAMAGED_RESIDUAL_STRAIN]);

		return eps_failure;
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
		if (rThisVariable == MP_EQUIVALENT_STRESS) rValue = mEquivalentStress;
		else if (rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN) rValue = mEquivalentPlasticStrain;
		else if (rThisVariable == MP_EQUIVALENT_PLASTIC_STRAIN_RATE) rValue = mEquivalentPlasticStrainRate;
		else if (rThisVariable == MP_DAMAGE) rValue = mDamage;
		else if (rThisVariable == MP_COMPACTION_RATIO) rValue = mAlpha;
		else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in RHT concrete 3D material law function GetValue double.";

		return(rValue);
	}

	void RHTConcrete3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in RHT concrete 3D material law function SetValue double.";
	}
} // Namespace Kratos