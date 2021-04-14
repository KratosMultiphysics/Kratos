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
		: JohnsonCookThermalPlastic3DLaw()
  { }

	JohnsonCookThermalPlastic2DPlaneStressLaw::JohnsonCookThermalPlastic2DPlaneStressLaw(const JohnsonCookThermalPlastic2DPlaneStressLaw& rOther)
  : JohnsonCookThermalPlastic3DLaw(rOther)
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
	  const Matrix identity = (GetStrainSize() == 3) ? IdentityMatrix(2) : IdentityMatrix(3);

	  // Convert vectors to matrices for easier manipulation
	  Matrix stress_old = Matrix(3, 3);
	  Matrix strain_increment = Matrix(3, 3);
	  MakeStrainStressMatrixFromVector((StrainVector - mStrainOld), strain_increment);
	  MakeStrainStressMatrixFromVector(StressVector, stress_old);

	  mStrainRate = std::sqrt(0.5 *
		  MPMStressPrincipalInvariantsUtility::CalculateMatrixDoubleContraction(strain_increment / CurrentProcessInfo[DELTA_TIME]));

	  // Material moduli
	  const double shear_modulus_G = (1.0 - mDamage) * MaterialProperties[YOUNG_MODULUS] / (2.0 + 2.0 * MaterialProperties[POISSON_RATIO]);
	  const double bulk_modulus_K = (1.0 - mDamage) * MaterialProperties[YOUNG_MODULUS] / (3.0 - 6.0 * MaterialProperties[POISSON_RATIO]);

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
	  Matrix stress_deviatoric_converged = Matrix(3, 3);
	  double delta_plastic_strain = 0.0;

	  // Plane stress - estimate through strain
	  double lame_lamda = (1.0 - mDamage)*MaterialProperties[YOUNG_MODULUS] * MaterialProperties[POISSON_RATIO]
		  / (1.0 + MaterialProperties[POISSON_RATIO]) / (1.0 - 2.0 * MaterialProperties[POISSON_RATIO]);
	  double delta_strain_33_old = lame_lamda * (strain_increment(0, 0) + strain_increment(1, 1))
		  / (lame_lamda + 3.0 * shear_modulus_G);
	  double pressure = stress_hydrostatic_new * -1.0;

	  if (j2_stress_trial > mYieldStressOld && mDamage < 0.95)
	  {
		  Matrix stress_trial = Matrix(3, 3);
		  stress_trial = stress_deviatoric_trial + stress_hydrostatic_new * identity;


		  double predicted_temperature = 0.0;
		  double yield_stress = 0.0;
		  double delta_pressure = 0.0;
		  double delta_hydrostatic = 0.0;
		  double sigma_33_old = stress_trial(2, 2);

		  const SizeType iteration_limit = 50;
		  const double tolerance = 1e-6;
		  bool is_converged = false;

		  double delta_strain_33_new = -1.0 * (strain_increment(0, 0) + strain_increment(1, 1));
		  if (std::abs(delta_strain_33_old - delta_strain_33_new) < 1e-9) 
			  delta_strain_33_new = 0.5 * strain_increment(1, 1);

		  while (!is_converged)
		  {
			  // eqn 16
			  delta_pressure = -1.0 * (bulk_modulus_K - 2.0 / 3.0 * shear_modulus_G)
				  * (strain_increment(0, 0) + strain_increment(1, 1) + delta_strain_33_new);
			  delta_hydrostatic = -1.0 * delta_pressure;

			  // eqn 17


		  }


		  delta_plastic_strain = (j2_stress_trial - mYieldStressOld) / (3.0 * shear_modulus_G + mHardeningMod);
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
	  if (delta_plastic_strain < 1e-12)
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

} // Namespace Kratos
