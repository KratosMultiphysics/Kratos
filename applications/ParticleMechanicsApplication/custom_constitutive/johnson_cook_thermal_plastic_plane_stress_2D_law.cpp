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

  void JohnsonCookThermalPlastic2DPlaneStressLaw::CalculateMaterialResponseKirchhoff(Kratos::ConstitutiveLaw::Parameters& rValues)
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
	  const Matrix identity = (GetStrainSize() == 3) ? IdentityMatrix(2) : IdentityMatrix(3);

	  // Convert vectors to matrices for easier manipulation
	  Matrix stress_old = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
	  Matrix strain_increment = (GetStrainSize() == 3) ? Matrix(2, 2) : Matrix(3, 3);
	  MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
	  MakeStrainStressMatrixFromVector(StressVector, stress_old);

	  mStrainRate = std::sqrt(0.5 *
		  MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(strain_increment / CurrentProcessInfo[DELTA_TIME]));

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
	  const double yield_stress_failure_ratio = 1e-3; // particle fails if current_yield/virgin_yield < yield_stress_failure_ratio
	  double delta_plastic_strain = 0.0;

	  if (j2_stress_trial > yield_stress&& mDamage < 0.95)
	  {
		  // Newton raphson setup
		  double gamma = mGammaOld;
		  double gamma_min = 0.0;
		  double gamma_max = j2_stress_trial / GetSqrt6() / shear_modulus_G;
		  bool is_converged = false;
		  const SizeType iteration_limit = 10000;
		  IndexType iteration = 0;
		  const double tolerance_delta_gamma = 1e-6;
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
			  // Particle has already failed. It can only take compressive volumetric stresses!
			  stress_deviatoric_converged.clear();
			  if (stress_hydrostatic_new > 0.0) stress_hydrostatic_new = 0.0;
		  }
		  else
		  {
			  stress_deviatoric_converged = stress_deviatoric_trial;
		  }
	  }

	  // Update equivalent stress
	  Matrix stress_converged = stress_deviatoric_converged + stress_hydrostatic_new * identity;
	  mEquivalentStress = std::sqrt(3.0 / 2.0 *
		  MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(stress_deviatoric_converged));

	  // Update damage
	  if (delta_plastic_strain > 1e-12)
	  {
		  // Evolution
		  if (mDamageInitiation < 1.0)
		  {
			  const double plastic_damage_onset = CalculateDamageOnsetPlasticStrain(stress_hydrostatic_new,
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
