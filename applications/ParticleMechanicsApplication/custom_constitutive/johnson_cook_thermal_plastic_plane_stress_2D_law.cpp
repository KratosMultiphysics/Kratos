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
#include "custom_constitutive/johnson_cook_thermal_plastic_plane_stress_2D_law.hpp"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
	JohnsonCookThermalPlastic2DPlaneStressLaw::JohnsonCookThermalPlastic2DPlaneStressLaw()
		: JohnsonCookThermalPlastic2DPlaneStrainLaw()
  { }

	JohnsonCookThermalPlastic2DPlaneStressLaw::JohnsonCookThermalPlastic2DPlaneStressLaw(const JohnsonCookThermalPlastic2DPlaneStressLaw& rOther)
  : JohnsonCookThermalPlastic2DPlaneStrainLaw(rOther)
  { }

  ConstitutiveLaw::Pointer JohnsonCookThermalPlastic2DPlaneStressLaw::Clone() const
  {
    return Kratos::make_shared<JohnsonCookThermalPlastic2DPlaneStressLaw>(*this);
  }

  JohnsonCookThermalPlastic2DPlaneStressLaw::~JohnsonCookThermalPlastic2DPlaneStressLaw()
  { }

  void JohnsonCookThermalPlastic2DPlaneStressLaw::CalculateMaterialResponseKirchhoffForwardEuler(Kratos::ConstitutiveLaw::Parameters& rValues)
  {
	  KRATOS_TRY

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
	  Vector& StressVector = rValues.GetStressVector();

	  // Convert vectors to matrices for easier manipulation
	  Matrix stress_old_2d = Matrix(2, 2);
	  Matrix strain_increment_2d = Matrix(2, 2);
	  MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment_2d);
	  MakeStrainStressMatrixFromVector(StressVector, stress_old_2d);

	  Matrix stress_old = ZeroMatrix(3);
	  Matrix strain_increment = ZeroMatrix(3);

	  for (size_t i = 0; i < 2; i++) {
		  for (size_t j = 0; j < 2; j++) {
			  stress_old(i, j) = stress_old_2d(i, j);
			  strain_increment(i, j) = strain_increment_2d(i, j);
		  }
	  }

	  mStrainRate = std::sqrt(0.5 *
		  MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(strain_increment / CurrentProcessInfo[DELTA_TIME]));

	  // Material moduli
	  const double shear_modulus_G = (1.0 - mDamage) * MaterialProperties[YOUNG_MODULUS] / (2.0 + 2.0 * MaterialProperties[POISSON_RATIO]);
	  const double bulk_modulus_K = (1.0 - mDamage) * MaterialProperties[YOUNG_MODULUS] / (3.0 - 6.0 * MaterialProperties[POISSON_RATIO]);
	  double lame_lamda = (1.0 - mDamage) * MaterialProperties[YOUNG_MODULUS] * MaterialProperties[POISSON_RATIO]
		  / (1.0 + MaterialProperties[POISSON_RATIO]) / (1.0 - 2.0 * MaterialProperties[POISSON_RATIO]);

	  // Follow equations from whir1988

	  // eq1
	  strain_increment(2, 2) = -1.0 * lame_lamda * (strain_increment(0, 0) + strain_increment(1, 1))
		  / (lame_lamda + 2.0 * shear_modulus_G);

	  // eq2
	  Matrix stress_trial = ZeroMatrix(3);
	  for (size_t i = 0; i < 3; ++i) {
		  for (size_t j = 0; j < 3; ++j) {
			  if (i != j) stress_trial(i, j) = stress_old(i,j) + 2.0 * shear_modulus_G * strain_increment(i, j);
		  }
	  }

	  // eq3
	  Matrix stress_xi_trial = stress_trial;

	  // eq 5
	  double delta_pressure = -1.0 * (bulk_modulus_K - 2.0 / 3.0 * shear_modulus_G)
		  * MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(strain_increment);

	  // eq 6
	  stress_trial(0, 0) = stress_old(0, 0) + 2.0 * shear_modulus_G * strain_increment(0, 0) - delta_pressure;
	  stress_trial(1, 1) = stress_old(1, 1) + 2.0 * shear_modulus_G * strain_increment(1, 1) - delta_pressure;
	  stress_trial(2, 2) = stress_old(2, 2) + 2.0 * shear_modulus_G * strain_increment(2, 2) - delta_pressure;

	  // eq7
	  double pressure = -1.0/3.0* MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_trial);

	  // eq8
	  stress_xi_trial(0, 0) = stress_trial(0, 0) + pressure;
	  stress_xi_trial(1, 1) = stress_trial(1, 1) + pressure;
	  // eq9
	  stress_xi_trial(2, 2) = -1.0* stress_xi_trial(0, 0) - stress_xi_trial(1, 1);

	  const double j2_stress_trial = std::sqrt(3.0 / 2.0 *
		  MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_xi_trial));

	  // Declare deviatoric stress matrix to be used later
	  Matrix stress_converged = Matrix(2, 2);
	  double delta_plastic_strain = 0.0;

	  if (j2_stress_trial > mYieldStressOld && mDamage < 0.95)
	  {
		  // eq12
		  stress_trial(2, 2) = 2.0 * shear_modulus_G * strain_increment(2, 2) - delta_pressure;

		  double predicted_temperature = 0.0;
		  double yield_stress = 0.0;
		  double delta_pressure = 0.0;
		  double delta_hydrostatic = 0.0;
		  double j2_stress = 0.0;
		  double sigma_33_old = stress_trial(2, 2);

		  const SizeType iteration_limit = 50;
		  SizeType iteration = 1;
		  const double tolerance = 1e-6;
		  bool is_converged = false;

		  // eqn 13
		  double delta_strain_33_i_minus_1 = strain_increment(2, 2);

		  // eqn 14
		  double delta_strain_33_i = -1.0 * (strain_increment(0, 0) + strain_increment(1, 1));
		  if (std::abs(delta_strain_33_i_minus_1 - delta_strain_33_i) < 1e-9)
			  delta_strain_33_i = 0.5 * strain_increment(1, 1);

		  double delta_strain_33_i_plus_1 = 0.0;

		  while (!is_converged)
		  {
			  // eqn 16
			  delta_pressure = -1.0 * (bulk_modulus_K - 2.0 / 3.0 * shear_modulus_G)
				  * (strain_increment(0, 0) + strain_increment(1, 1) + delta_strain_33_i);

			  // eqn 17
			  stress_trial(0, 0) = stress_old(0, 0) + 2.0 * shear_modulus_G * strain_increment(0, 0) - delta_pressure;
			  stress_trial(1, 1) = stress_old(1, 1) + 2.0 * shear_modulus_G * strain_increment(1, 1) - delta_pressure;

			  // eqn18
			  pressure = -1.0 / 3.0 * MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_trial);

			  // eqn19
			  stress_xi_trial(0, 0) = stress_trial(0, 0) + pressure;
			  stress_xi_trial(1, 1) = stress_trial(1, 1) + pressure;

			  //eqn20
			  stress_xi_trial(2, 2) = -1.0 * stress_xi_trial(0, 0) - stress_xi_trial(1, 1);

			  //eqn21
			  j2_stress = std::sqrt(3.0 / 2.0 *
				  MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_xi_trial));
			  delta_plastic_strain = (j2_stress - mYieldStressOld) / (3.0 * shear_modulus_G + mHardeningMod);

			  // eqn22
			  sigma_33_old = stress_trial(2, 2);
			  stress_trial(2, 2) = sigma_33_old - 3.0 * shear_modulus_G * delta_plastic_strain * stress_xi_trial(2, 2) / j2_stress;

			  // eqn23
			  delta_strain_33_i_plus_1 = delta_strain_33_i_minus_1 - sigma_33_old * (delta_strain_33_i - delta_strain_33_i_minus_1)
				  / (1.0e-6 + stress_trial(2, 2) - sigma_33_old);

			  // eqn24
			  double residual = std::abs(delta_strain_33_i - delta_strain_33_i_minus_1) / std::abs(delta_strain_33_i_plus_1);
			  if (residual < 1e-6)
			  {
				  is_converged = true;

				  strain_increment(2, 2) = delta_strain_33_i_plus_1;

				  stress_converged = stress_trial - 3.0 * shear_modulus_G * delta_plastic_strain * stress_xi_trial / j2_stress;
				  stress_converged(2, 2) = 0.0;
				  stress_converged.resize(2, 2, true);

				  // Normal plasticity updating
				  mEquivalentStress = j2_stress;
				  mEquivalentPlasticStrainOld += delta_plastic_strain;
				  mPlasticStrainRateOld = delta_plastic_strain / CurrentProcessInfo[DELTA_TIME];

				  // Use only old yield stress for temperature factor in new yield stress
				  predicted_temperature = mTemperatureOld + MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT]
					  / MaterialProperties[DENSITY] / MaterialProperties[SPECIFIC_HEAT]
					  * (mYieldStressOld)*delta_plastic_strain;
				  yield_stress = CalculateHardenedYieldStress(MaterialProperties,
					  mEquivalentPlasticStrainOld, mPlasticStrainRateOld, predicted_temperature);

				  // Now update the temperature with the new yield stress
				  mTemperatureOld += MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] / MaterialProperties[DENSITY]
					  / MaterialProperties[SPECIFIC_HEAT] * (yield_stress)*delta_plastic_strain;

				  // Finalize plasticity
				  mHardeningMod = (yield_stress - mYieldStressOld) / delta_plastic_strain; // already includes damage
				  mEnergyDissipated += (yield_stress) / delta_plastic_strain / MaterialProperties[DENSITY];
				  mYieldStressOld = yield_stress;

				  break;
			  }
			  else
			  {
				  iteration += 1;
				  if (iteration > iteration_limit)
				  {
					  KRATOS_INFO("ADF") << "INTERATION LIMIT EXCEEDED!";
					  KRATOS_ERROR << "ERROR";
				  }
				  else
				  {
					  delta_strain_33_i_minus_1 = delta_strain_33_i;
					  delta_strain_33_i = delta_strain_33_i_plus_1;
				  }
			  }
		  }
	  }
	  else
	  {
			mEquivalentStress = j2_stress_trial;

			// eqn 11
			for (size_t i = 0; i < 2; i++)
			{
				for (size_t j = 0; j < 2; j++)
				{
					stress_converged(i, j) = stress_old(i, j) + 2.0 * shear_modulus_G * strain_increment(i, j);
					if (i == j) stress_converged(i, j) -= delta_pressure;
				}
			}

		  if (mDamage > 0.95)
		  {
			  stress_converged.clear();
			  // Particle has already failed. It can only take compressive volumetric stresses!
			  if (pressure > 0.0) {
				  // compressive stress
				  stress_converged(0, 0) -= pressure;
				  stress_converged(1, 1) -= pressure;
			  }
		  }

	  }

	  // Update damage
	  if (delta_plastic_strain < 1e-12)
	  {
		  // Evolution
		  if (mDamageInitiation < 1.0)
		  {
			  const double plastic_damage_onset = CalculateDamageOnsetPlasticStrain(-1.0*pressure,
				  mEquivalentStress, mPlasticStrainRateOld, mTemperatureOld, MaterialProperties);

			  mDamageInitiation += (delta_plastic_strain / plastic_damage_onset);
		  }
		  else
		  {
			  // True material damage
			  if (mDamageInitiationDisp < 1e-12)
			  {
				  // We have just started damage - store the current plastic displacement
				  mDamageInitiationDisp = mCharLength * mEquivalentPlasticStrainOld;
				  mFailureDisp = 2.0 * mFractureEnergy / mYieldStressOld; //value of the yield stress at the time when the failure criterion is reached
				  mDamage = 1e-9;
				  if (mFractureEnergy < 1e-9) mDamage = 1.0;
			  }

			  if (mDamage < 1.0)
			  {
				  // Softening regularised with Hillerborg approach
				  const bool is_exponential_softening = true;
				  const double plastic_disp_after_onset = mCharLength * mEquivalentPlasticStrainOld - mDamageInitiationDisp;
				  if (is_exponential_softening)
				  {
					  // https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-damageevolductile.htm
					  const double undamaged_yield = mYieldStressOld / (1.0 - mDamage);
					  mDamage = 1.0 - std::exp(-1.0 * plastic_disp_after_onset * undamaged_yield / mFractureEnergy);
				  }
				  else
				  {
					  // simple linear damage evolution
					  if (mFailureDisp - mDamageInitiationDisp < 1e-9)  mDamage = 1.0;
					  else mDamage = plastic_disp_after_onset / (mFailureDisp - mDamageInitiationDisp); // simple linear damage law
				  }
				  mDamage = std::min(1.0, mDamage);
				  mDamage = std::max(0.0, mDamage);
				  if (mDamage > 0.95) mDamage = 1.0;
			  }
		  }
	  }


	  // Store stresses and strains
	  MakeStrainStressVectorFromMatrix(stress_converged, StressVector);
	  mStrainOld = StrainVector;

	  // Udpdate internal energy
	  for (size_t i = 0; i < stress_converged.size1(); ++i) {
		  for (size_t j = 0; j < stress_converged.size2(); ++j) {
			  mEnergyInternal += 0.5 * (stress_converged(i, j) + stress_old(i, j)) * strain_increment(i, j) / MaterialProperties[DENSITY];
		  }
	  }
	  KRATOS_CATCH("")
  }

  void JohnsonCookThermalPlastic2DPlaneStressLaw::CalculateMaterialResponseKirchhoffBackwardEuler(Kratos::ConstitutiveLaw::Parameters& rValues)
  {
	  KRATOS_TRY

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
	  Vector& StressVector = rValues.GetStressVector();
	  const Matrix identity = IdentityMatrix(2);

	  // Convert vectors to matrices for easier manipulation
	  Matrix stress_old = Matrix(2, 2);
	  Matrix strain_increment = Matrix(2, 2);
	  MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
	  MakeStrainStressMatrixFromVector(StressVector, stress_old);

	  mStrainRate = std::sqrt(0.5 *
		  MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(strain_increment / CurrentProcessInfo[DELTA_TIME]));

	  // Material moduli
	  const double shear_modulus_G = (1.0 - mDamage) * MaterialProperties[YOUNG_MODULUS] / (2.0 + 2.0 * MaterialProperties[POISSON_RATIO]);
	  const double bulk_modulus_K = (1.0 - mDamage) * MaterialProperties[YOUNG_MODULUS] / (3.0 - 6.0 * MaterialProperties[POISSON_RATIO]);

	  // Predict plane stresses
	  Matrix stress_converged = stress_old;
	  double temp = (1.0-mDamage)* MaterialProperties[YOUNG_MODULUS] / (1.0 - MaterialProperties[POISSON_RATIO]* MaterialProperties[POISSON_RATIO]);
	  stress_converged += temp * (1.0 - MaterialProperties[POISSON_RATIO]) * strain_increment;
	  stress_converged += temp * MaterialProperties[POISSON_RATIO] * identity * (strain_increment(0, 0) + strain_increment(1, 1));

	  // Calc j2 stress
	  double j2_stress_trial = stress_converged(0, 0) * stress_converged(0, 0)
		  + stress_converged(1, 1) * stress_converged(1, 1)
		  - stress_converged(0, 0) * stress_converged(1, 1)
		  + 3.0 * stress_converged(0, 1) * stress_converged(0, 1);
	  j2_stress_trial = std::sqrt(j2_stress_trial);

	  double stress_hydrostatic_new = MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_converged) / 3.0;

	  // Declare deviatoric stress matrix to be used later

	  // Assume current yield stress is the same as the old (elastic predictor)
	  double yield_stress = mYieldStressOld;
	  const double yield_stress_failure_ratio = 1e-3; // particle fails if current_yield/virgin_yield < yield_stress_failure_ratio
	  double delta_plastic_strain = 0.0;

	  if (j2_stress_trial > yield_stress && mDamage < 0.95)
	  {
		  // Recast as axisymmetric problem and secant iteration for epsilon_33 [plane stress plasticity]
		  double epsilon_33_n_minus_1 = (StrainVector[0] + StrainVector[1]) * -1.0 * (bulk_modulus_K - 2.0 / 3.0 * shear_modulus_G) / (4.0 / 3.0 * shear_modulus_G + bulk_modulus_K); // [plane_stress_secant]
		  double epsilon_33_n = (StrainVector[0] + StrainVector[1]) * -1.0 * MaterialProperties[POISSON_RATIO] / (1.0 - MaterialProperties[POISSON_RATIO]); //[plane_stress_secant]
		  double epsilon_33_n_plus_1 = 0.0;

		  double stress_33_n_minus_1, stress_33_n, stress_33_n_plus_1, epsilon_33;

		  bool plane_stress_converged = false;
		  const IndexType plane_stress_iteration_limit = 11;
		  IndexType plane_stress_iteration = 1;
		  Matrix strain_increment_axisym = ZeroMatrix(3, 3);
		  const Matrix identity_axisym = IdentityMatrix(3);
		  Matrix stress_old_axisym = ZeroMatrix(3, 3);
		  Matrix stress_deviatoric_converged;

		  double yield_function, yield_function_gradient, dYield_dGamma,
			  delta_gamma, predicted_eps, predicted_eps_rate,
			  predicted_temperature, gamma;

		  while (!plane_stress_converged)
		  {
			  // Set axisym
			  if (plane_stress_iteration == 1) epsilon_33 = epsilon_33_n_minus_1;
			  else epsilon_33 = epsilon_33_n;
			  for (size_t i = 0; i < 2; i++) {
				  for (size_t j = 0; j < 2; j++) {
					  strain_increment_axisym(i, j) = strain_increment(i, j);
					  stress_old_axisym(i, j) = stress_old(i, j);
				  }
			  }
			  strain_increment_axisym(2, 2) = epsilon_33;

			  Matrix strain_increment_deviatoric = strain_increment_axisym -MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(strain_increment_axisym) / 3.0 * identity_axisym;
			  Matrix stress_deviatoric_old = stress_old_axisym -  MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_old_axisym) / 3.0 * identity_axisym;

			  // Calculate trial (predicted) j2 stress
			  stress_hydrostatic_new = MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(stress_old_axisym) / 3.0 + bulk_modulus_K * MPMStressPrincipalInvariantsUtility::CalculateMatrixTrace(strain_increment_axisym);
			  Matrix stress_deviatoric_trial = stress_deviatoric_old + 2.0 * shear_modulus_G * strain_increment_deviatoric;
			  j2_stress_trial = std::sqrt(3.0 / 2.0 *MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_deviatoric_trial));
			  stress_deviatoric_converged = stress_deviatoric_trial;

			  if (j2_stress_trial > yield_stress)
			  {
				  // Newton raphson setup
				  gamma = mGammaOld;
				  double gamma_min = 0.0;
				  double gamma_max = j2_stress_trial / GetSqrt6() / shear_modulus_G;
				  bool is_converged = false;
				  const SizeType iteration_limit = 10000;
				  IndexType iteration = 0;
				  const double tolerance_delta_gamma = 1e-6;

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
						  MaterialProperties, predicted_eps, predicted_eps_rate, predicted_temperature) / CurrentProcessInfo[DELTA_TIME];
					  dYield_dGamma += MaterialProperties[TAYLOR_QUINNEY_COEFFICIENT] * yield_stress / MaterialProperties[DENSITY] /
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
						  if (gamma < 1e-6) is_converged = true;
						  else
						  {
							#pragma omp critical
							  {
								  KRATOS_INFO("Johnson Cook Material Model") << " Johnson Cook iteration limit exceeded"
									  << "\ngamma = " << gamma
									  << "\ndamage = " << mDamage
									  << "\ndamage initiation = " << mDamageInitiation
									  << "\nplastic strain old = " << mEquivalentPlasticStrainOld
									  << "\nyield stress / virgin yield = " << yield_stress / mYieldStressVirgin
									  << "\ndelta_gamma = " << delta_gamma
									  << "\nyield_function" << yield_function
									  << "\n\n\n";

								  KRATOS_ERROR << "Johnson Cook iteration limit exceeded";
							  }
						  }
					  }
				  }
				  // Correct trial stress
				  stress_deviatoric_converged = stress_deviatoric_trial - 2.0 * shear_modulus_G * gamma * flow_direction_normalized;
			  }


			  // Plane stress check -------------------------------------------------
			  Matrix total_stress = stress_deviatoric_converged + stress_hydrostatic_new * identity_axisym;

			  if (plane_stress_iteration == 1) stress_33_n_minus_1 = total_stress(2, 2);

			  if (plane_stress_iteration > 1)
			  {
				  stress_33_n = total_stress(2, 2);
				  epsilon_33_n_plus_1 = epsilon_33_n - stress_33_n * (epsilon_33_n - epsilon_33_n_minus_1) / (stress_33_n - stress_33_n_minus_1);

				  double plane_stress_residual = std::abs(epsilon_33_n - epsilon_33_n_minus_1) / std::abs(epsilon_33_n_plus_1);
				  if (plane_stress_residual < 1e-4)
				  {
					  // converged
					  stress_converged = total_stress;
					  stress_converged.resize(2, 2, true);
					  mEquivalentStress = j2_stress_trial;
					  plane_stress_converged = true;
					  break;
				  }
				  else
				  {
					  if (plane_stress_iteration > plane_stress_iteration_limit)
					  {
						#pragma omp critical
						  {
							  KRATOS_INFO("Johnson Cook Material Model") << " Plane stress iteration limit exceeded"
								  << "\nstress_33_n_minus_1 = " << stress_33_n_minus_1
								  << "\nstress_33_n = " << stress_33_n
								  << "\n\n\n";

							  KRATOS_ERROR << "Johnson Cook iteration limit exceeded";
						  }
					  }
					  else
					  {
						  // move all data one step older
						  epsilon_33_n_minus_1 = epsilon_33_n;
						  epsilon_33_n = epsilon_33_n_plus_1;
						  stress_33_n_minus_1 = stress_33_n;


					  }


				  }
			  }
			  plane_stress_iteration += 1;
		  }


		  // Store plastic deformation quantities
		  delta_plastic_strain = predicted_eps - mEquivalentPlasticStrainOld;
		  mEnergyDissipated += gamma / GetSqrt6() / MaterialProperties[DENSITY] * (mYieldStressOld + yield_stress);
		  mEquivalentPlasticStrainOld = predicted_eps;
		  mPlasticStrainRateOld = predicted_eps_rate;
		  mTemperatureOld = predicted_temperature;
		  mYieldStressOld = yield_stress;
	  }
	  else
	  {
		  if (mDamage > 0.95)
		  {
			  // Particle has already failed.
			  stress_converged.clear();
			  mEquivalentStress = 0.0;
		  }
		  else
		  {
			  mEquivalentStress = j2_stress_trial;
		  }
	  }

	  // Update damage
	  if (delta_plastic_strain > 1e-12)
	  {
		  // Evolution
		  if (mDamageInitiation < 1.0)
		  {
			  const double plastic_damage_onset = CalculateDamageOnsetPlasticStrain(stress_hydrostatic_new,
				  mEquivalentStress, mPlasticStrainRateOld, mTemperatureOld, MaterialProperties);

			  if (plastic_damage_onset < 1e-12) mDamageInitiation = 1.0;
			  else mDamageInitiation += (delta_plastic_strain / plastic_damage_onset);
		  }
		  else
		  {
			  // True material damage
			  if (mDamageInitiationDisp < 1e-12)
			  {
				  // We have just started damage - store the current plastic displacement
				  mDamageInitiationDisp = mCharLength * mEquivalentPlasticStrainOld;
				  mFailureDisp = 2.0 * mFractureEnergy / yield_stress; //value of the yield stress at the time when the failure criterion is reached
				  mDamage = 1e-9;
				  if (mFractureEnergy < 1e-9) mDamage = 1.0;
			  }

			  if (mDamage < 1.0)
			  {
				  // Softening regularised with Hillerborg approach
				  const bool is_exponential_softening = true;
				  const double plastic_disp_after_onset = mCharLength * mEquivalentPlasticStrainOld - mDamageInitiationDisp;
				  if (is_exponential_softening)
				  {
					  // https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-damageevolductile.htm
					  const double undamaged_yield = yield_stress / (1.0 - mDamage);
					  mDamage = 1.0 - std::exp(-1.0 * plastic_disp_after_onset * undamaged_yield / mFractureEnergy);
				  }
				  else
				  {
					  // simple linear damage evolution
					  if (mFailureDisp - mDamageInitiationDisp < 1e-9)  mDamage = 1.0;
					  else mDamage = plastic_disp_after_onset / (mFailureDisp - mDamageInitiationDisp); // simple linear damage law
				  }
				  mDamage = std::min(1.0, mDamage);
				  mDamage = std::max(0.0, mDamage);
				  if (mDamage > 0.95) mDamage = 1.0;
			  }
		  }
	  }


	  // Store stresses and strains
	  MakeStrainStressVectorFromMatrix(stress_converged, StressVector);
	  mStrainOld = StrainVector;

	  // Udpdate internal energy
	  for (size_t i = 0; i < stress_converged.size1(); ++i) {
		  for (size_t j = 0; j < stress_converged.size2(); ++j) {
			  mEnergyInternal += 0.5 * (stress_converged(i, j) + stress_old(i, j)) * strain_increment(i, j) / MaterialProperties[DENSITY];
		  }
	  }
	  KRATOS_CATCH("")
  }

} // Namespace Kratos
