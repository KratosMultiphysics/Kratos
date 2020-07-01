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

		mStrainOld = ZeroVector(GetStrainSize());
		mDamage = 0.0;
		mAlpha = 1.1;

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
		const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
		CheckIsExplicitTimeIntegration(CurrentProcessInfo);
		const Properties& MaterialProperties = rValues.GetMaterialProperties();


		// Get old stress vector and current strain vector
		const Vector StrainVector = rValues.GetStrainVector();
		Vector& StressVector = rValues.GetStressVector();
		const Matrix identity = (GetStrainSize() == 3) ? IdentityMatrix(2) : IdentityMatrix(3);

		// Convert vectors to matrices for easier manipulation
		Matrix stress_old = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		Matrix strain_increment = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		Matrix strain_new = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
		MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
		MakeStrainStressMatrixFromVector(StressVector, stress_old);
		MakeStrainStressMatrixFromVector(StrainVector, strain_new);

		// Calculate deviatoric quantities
		const Matrix strain_increment_deviatoric = strain_increment -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(strain_increment) / 3.0 * identity;
		const Matrix stress_deviatoric_old = stress_old -
			MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_old) / 3.0 * identity;
		Matrix stress_deviatoric_trial = stress_deviatoric_old + 2.0 * MaterialProperties[SHEAR_MODULUS] * strain_increment_deviatoric;

		// Calculate pressure from Equation of State
		const double pressure = CalculatePressureFromEOS(MaterialProperties, strain_new, stress_old);

		// Prepare elastic predictor
		const double eff_stress_trial = std::sqrt(3.0 / 2.0 *
			MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_deviatoric_trial));
		const double lode_angle =
			MPMStressPrincipalInvariantsUtility::CalculateLodeAngleFromDeviatoricStressTensor(stress_deviatoric_trial);

		double limit_elastic = CalculateElasticLimit(const double lode_angle, const double eps_rate, const double eps, const Properties & MaterialProperties)



		const double density_current = (MaterialProperties [DENSITY]) / rValues.GetDeterminantF();


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
		bool alpha_smoothing = false; // used to prevent stick-slipping in solving alpha
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
				rMaterialProperties[RHT_EOS_ALPHA0] / rMaterialProperties[REFERENCE_DENSITY] - 1.0;

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

			if (std::abs(pressure-pressure_old) < rMaterialProperties[RHT_COMPRESSIVE_STRENGTH]/1e6)
			{
				is_converged = true;
				break;
			}
			else
			{
				iteration += 1;
				if (iteration > iteration_limit)
				{
					KRATOS_INFO("CalculatePressureFromEOS") << "Iteration = " << iteration << "\n"
						<< "Pressure = " << pressure << "\n"
						<< "Pressure_old = " << pressure_old << "\n"
						<< "Alpha_trial = " << alpha_trial << "\n"
						<< "Alpha_old = " << alpha_old << "\n";
					KRATOS_ERROR << "CalculatePressureFromEOS iteration limit exceeded\n";
				}
				if (iteration > 5) alpha_smoothing = true;
				pressure_old = pressure;
				alpha_old = alpha_trial;
			}
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
		const double EPSrate, const double EPS, const Properties& rMaterialProperties)
	{
		// ref[1] eqn 14
		const double pressure_star = Pressure / rMaterialProperties[RHT_COMPRESSIVE_STRENGTH];
		const double pore_crush_star = mPoreCrushPressure / rMaterialProperties[RHT_COMPRESSIVE_STRENGTH];

		const double fe_elastic_factor = CalculateFeElasticFactor(pressure_star,
			EPSrate, rMaterialProperties);

		const double fc_cap_factor = CalculateFcCapFactor(pressure_star, pore_crush_star,
			EPSrate, EPS, rMaterialProperties);

		const double fr_rate_factor = CalculateFrRateFactor(pressure_star, EPSrate,
			rMaterialProperties);

		const double normalized_yield = CalculateNormalizedYield(pressure_star/ fe_elastic_factor,
			fr_rate_factor, EPSrate, rMaterialProperties);

		const double q_deviatoric_shape_factor = CalculateDeviatoricShapeFactorQ(pressure_star, rMaterialProperties);
		const double r_triaxiality = CalculateTriaxialityR(LodeAngle, q_deviatoric_shape_factor);

		return rMaterialProperties[RHT_COMPRESSIVE_STRENGTH] * normalized_yield * r_triaxiality *
			fe_elastic_factor * fc_cap_factor;
	}

	const double RHTConcrete3DLaw::CalculateFeElasticFactor(const double PressureStar,
		const double EPSrate, const Properties& rMaterialProperties)
	{
		// ref[1] eqn15
		const double gc_star = rMaterialProperties[RHT_GC_STAR];
		const double gt_star = rMaterialProperties[RHT_GT_STAR];
		const double ft_star = rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH];

		double fe_elastic_factor = 0.0;
		CalculateRateFactors(EPSrate, rMaterialProperties);

		if (3.0* PressureStar >= mRateFactorCompression* gc_star)
		{
			fe_elastic_factor = gc_star;
		}
		else if (3.0 * PressureStar >= -1.0* mRateFactorTension* gt_star * ft_star)
		{
			fe_elastic_factor = gc_star - (gt_star - gc_star) * (3.0 * PressureStar - mRateFactorCompression * gc_star)/
				(mRateFactorCompression * gc_star + mRateFactorTension * gt_star * ft_star);
		}
		else
		{
			fe_elastic_factor = gt_star;
		}

		return fe_elastic_factor;
	}


	const double RHTConcrete3DLaw::CalculateFcCapFactor(const double PressureStar,
		const double PoreCrushStar, const double EPSrate, const double EPS,
		const Properties& rMaterialProperties)
	{
		//ref [1] eqn16
		double fc_cap_factor;
		const double degraded_shear_mod = rMaterialProperties[SHEAR_MODULUS] *
			rMaterialProperties[RHT_SHEAR_MOD_REDUCTION_FACTOR];
		CalculateRateFactors(EPSrate, rMaterialProperties);
		const double initial_pressure_star =
			mRateFactorCompression * rMaterialProperties[RHT_GC_STAR] / 3.0 +
			degraded_shear_mod * EPS / rMaterialProperties[RHT_COMPRESSIVE_STRENGTH];

		if (PressureStar >= PoreCrushStar) fc_cap_factor = 0.0;
		else if (std::abs(PressureStar - PoreCrushStar) < 1e-6) fc_cap_factor = 0.0;
		else if (PressureStar >= initial_pressure_star)
		{
			fc_cap_factor = std::sqrt(1.0 - std::pow((PressureStar - initial_pressure_star) /
			(PoreCrushStar - initial_pressure_star), 2.0));
		}
		else fc_cap_factor = 1.0;

		return fc_cap_factor;
	}


	const double RHTConcrete3DLaw::CalculateFrRateFactor(const double PressureStarForRate,
		const double EPSRate, const Properties& rMaterialProperties)
	{
		// ref[1] eqn10
		CalculateRateFactors(EPSRate, rMaterialProperties);
		double fr_rate_factor;
		if (3.0 * PressureStarForRate >= mRateFactorCompression) fr_rate_factor = mRateFactorCompression;
		else if (3.0 * PressureStarForRate >= -1.0 * mRateFactorTension * rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH])
		{
			fr_rate_factor = mRateFactorCompression - (mRateFactorTension - mRateFactorCompression) *
				(3.0 * PressureStarForRate - mRateFactorCompression) /
				(mRateFactorCompression + mRateFactorTension * rMaterialProperties[RHT_RELATIVE_TENSILE_STRENGTH]);
		}
		else fr_rate_factor = mRateFactorTension;

		return fr_rate_factor;
	}


	const double RHTConcrete3DLaw::CalculateNormalizedYield(const double PressureStar,
		const double RateFactor, const double EPSrate, const Properties& rMaterialProperties)
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


	void RHTConcrete3DLaw::CalculateRateFactors(const double EPSrate,
		const Properties& rMaterialProperties)
	{
		// ref[1] eqn11
		if (EPSrate == 0.0) {
			mRateFactorTension = 1.0;
			mRateFactorCompression = 1.0;
		}
		else
		{
			const double beta_t = 1.0 / (10.0 + 0.5 * rMaterialProperties[RHT_COMPRESSIVE_STRENGTH] / 1e6);
			const double beta_c = 1.0 / (5.0 + 0.75 * rMaterialProperties[RHT_COMPRESSIVE_STRENGTH] / 1e6);

			mRateFactorTension = std::pow(EPSrate/ rMaterialProperties[REFERENCE_TENSION_STRAIN_RATE], beta_t);
			mRateFactorCompression = std::pow(EPSrate/ rMaterialProperties[REFERENCE_COMPRESSION_STRAIN_RATE], beta_c);
		}
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