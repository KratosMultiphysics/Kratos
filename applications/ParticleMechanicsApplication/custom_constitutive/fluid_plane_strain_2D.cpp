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

#include "custom_constitutive/fluid_plane_strain_2D.hpp"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
	FluidPlaneStrain2DLaw::FluidPlaneStrain2DLaw()
		: HyperElastic3DLaw()
	{
	}


	FluidPlaneStrain2DLaw::FluidPlaneStrain2DLaw(const FluidPlaneStrain2DLaw& rOther)
		: HyperElastic3DLaw(rOther)
	{
	}


	ConstitutiveLaw::Pointer FluidPlaneStrain2DLaw::Clone() const
	{
		return Kratos::make_shared<FluidPlaneStrain2DLaw>(*this);
	}


	FluidPlaneStrain2DLaw::~FluidPlaneStrain2DLaw()
	{  }



	void FluidPlaneStrain2DLaw::InitializeMaterial(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const Vector& rShapeFunctionsValues)
	{
		KRATOS_TRY

		BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

		mStrainOld = ZeroVector(GetStrainSize());

		// Add static pressure
		const double rho = rMaterialProperties[DENSITY];
		double elevation = 0.0;
		const double integration_points_number = rElementGeometry.IntegrationPointsNumber();
		double weight = 1.0;
		const Matrix shape_functions = rElementGeometry.ShapeFunctionsValues();
		double check = 0.0;
		for (size_t int_p = 0; int_p < integration_points_number; int_p++)
		{
			if (integration_points_number > 1) weight = rElementGeometry.IntegrationPoints()[int_p].Weight();
			for (size_t i = 0; i < rElementGeometry.PointsNumber(); i++)
			{
				if (shape_functions(int_p, i) > 1e-12) elevation += weight * shape_functions(int_p, i) * rElementGeometry.GetPoint(i)[1];
			}
		}

		mPressure = 9.81*rho*(5.0-elevation);

		KRATOS_CATCH("")
	}


	void FluidPlaneStrain2DLaw::CalculateMaterialResponseKirchhoff(Kratos::ConstitutiveLaw::Parameters& rValues)
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
		const double detF = rValues.GetDeterminantF();

		// Convert vectors to matrices for easier manipulation
		const Vector strain_increment_rate_vec = (StrainVector - mStrainOld) / CurrentProcessInfo[DELTA_TIME];
		Matrix strain_increment_rate_mat = (GetStrainSize() == 3) ? ZeroMatrix(2, 2) : ZeroMatrix(3, 3);
			// Normal components
		strain_increment_rate_mat(0, 0) = strain_increment_rate_vec[0];
		strain_increment_rate_mat(1, 1) = strain_increment_rate_vec[1];
			// Shear components
		strain_increment_rate_mat(0, 1) = 0.5 * strain_increment_rate_vec[2]; //xy
		strain_increment_rate_mat(1, 0) = strain_increment_rate_mat(0, 1);

		// Material properties
		const double bulk_modulus_K = MaterialProperties[BULK_MODULUS];
		const double dynamic_viscosity = MaterialProperties[VISCOSITY];
		mPressure = bulk_modulus_K * (1.0 - detF) / detF;

		Matrix stress_mat = -1.0* mPressure * identity + 2.0 * dynamic_viscosity * strain_increment_rate_mat; // [mast2013]

		// Store stresses and strains
		StressVector.clear();
		StressVector[0] = stress_mat(0, 0);
		StressVector[1] = stress_mat(1, 1);
		StressVector[2] = 2.0 * stress_mat(0, 1); //xy
		mStrainOld = StrainVector;

		KRATOS_CATCH("")
	}


	int FluidPlaneStrain2DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

		KRATOS_ERROR_IF(BULK_MODULUS.Key()==0 || rMaterialProperties[BULK_MODULUS] < 0.0) << "BULK_MODULUS has key zero or invalid value" << std::endl;
		KRATOS_ERROR_IF(DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.00) << "DENSITY has Key zero or invalid value " << std::endl;
		KRATOS_ERROR_IF(VISCOSITY.Key() == 0 || rMaterialProperties[VISCOSITY] < 0.00) << "VISCOSITY has Key zero or invalid value " << std::endl;

		//if (check_base > 1) return 1;
		return 0;

		KRATOS_CATCH("")
	}


	bool FluidPlaneStrain2DLaw::CheckParameters(Parameters& rValues)
	{
		return rValues.CheckAllParameters();
	}


	void FluidPlaneStrain2DLaw::CheckIsExplicitTimeIntegration(const ProcessInfo& rCurrentProcessInfo)
	{
		const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
			? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
			: false;
		KRATOS_ERROR_IF_NOT(is_explicit) << "The Johnson Cook MPM material law is currently limited to explicit time integration only";
	}



	void FluidPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
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


	double& FluidPlaneStrain2DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
	{
		if (rThisVariable == PRESSURE)
		{
			rValue = mPressure;
		}
		else KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in FluidPlaneStrain2DLaw material law function GetValue double.";

		return(rValue);
	}


	void FluidPlaneStrain2DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR << "Variable " << rThisVariable << " not implemented in FluidPlaneStrain2DLaw material law function SetValue double.";
	}
} // Namespace Kratos