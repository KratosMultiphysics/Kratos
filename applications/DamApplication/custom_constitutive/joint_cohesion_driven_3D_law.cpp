//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Javier San Mauro Saiz
//                   Joaquin Irazabal Gonzalez
//

// Application includes
#include "custom_constitutive/joint_cohesion_driven_3D_law.hpp"

namespace Kratos
{

int JointCohesionDriven3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    // Verify ProcessInfo variables
    KRATOS_CHECK_VARIABLE_KEY(IS_CONVERGED);

    // Verify Properties variables
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    if(rMaterialProperties.Has(YOUNG_MODULUS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YOUNG_MODULUS not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO);
    if(rMaterialProperties.Has(POISSON_RATIO)) {
        KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] < -1.0) << "POISSON_RATIO has an invalid value lower than -1.0" << std::endl;
        KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] >= 0.5) << "POISSON_RATIO has an invalid value greater or equal to 0.5 " << std::endl;
    } else {
        KRATOS_ERROR << "POISSON_RATIO not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(FRICTION_COEFFICIENT);
    if(rMaterialProperties.Has(FRICTION_COEFFICIENT)) {
        KRATOS_ERROR_IF(rMaterialProperties[FRICTION_COEFFICIENT] < 0.0) << "FRICTION_COEFFICIENT has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "FRICTION_COEFFICIENT not defined" << std::endl;
    }

        KRATOS_CHECK_VARIABLE_KEY(COHESION);
    if(rMaterialProperties.Has(COHESION)) {
        KRATOS_ERROR_IF(rMaterialProperties[COHESION] < 0.0) << "COHESION has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "COHESION not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------

void JointCohesionDriven3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable = 1.0;
}

//----------------------------------------------------------------------------------------

void JointCohesionDriven3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        rValues.CheckAllParameters();

        ConstitutiveLawVariables Variables;
        this->InitializeConstitutiveLawVariables(Variables,rValues);

        this->ComputeEquivalentStrain(Variables,rValues);
        this->CheckLoadingFunction(Variables,rValues);

        if(Variables.LoadingFlag)
        {
            mStateVariable = Variables.EquivalentStrain;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void JointCohesionDriven3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                 Parameters& rValues)

{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    rVariables.YoungModulus = MaterialProperties[YOUNG_MODULUS];
    rVariables.PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    rVariables.Cohesion = MaterialProperties[COHESION];
    rVariables.YieldStress = rVariables.YoungModulus;
    rVariables.FrictionCoefficient = MaterialProperties[FRICTION_COEFFICIENT];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void JointCohesionDriven3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();
	//Triangular area
    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        rVariables.EquivalentStrain = 1.0;
        if (mStateVariable == 1.0)
        {
            double tau0 = rVariables.YieldStress * StrainVector[0];
            double tau1 = rVariables.YieldStress * StrainVector[1];
		    double sigma = rVariables.YoungModulus * StrainVector[2];

		    double broken_limit = (-rVariables.FrictionCoefficient * rVariables.YoungModulus * StrainVector[2])+ rVariables.Cohesion;

		    if (sigma > (rVariables.Cohesion/rVariables.FrictionCoefficient)) rVariables.EquivalentStrain = 0.0;
		    if (abs (tau0) > broken_limit) rVariables.EquivalentStrain = 0.0;
		    if (abs (tau1) > broken_limit) rVariables.EquivalentStrain = 0.0;
        }
    }
    else // Contact between interfaces
    {
        rVariables.EquivalentStrain = 1.0;
        if (mStateVariable == 1.0)
        {
		    double tau0 = rVariables.YieldStress * StrainVector[0];
            double tau1 = rVariables.YieldStress * StrainVector[1];
		    double sigma = rVariables.YoungModulus * StrainVector[2];

		    double broken_limit = (-rVariables.FrictionCoefficient * rVariables.YoungModulus * StrainVector[2])+ rVariables.Cohesion;

		    if (sigma > (rVariables.Cohesion/rVariables.FrictionCoefficient)) rVariables.EquivalentStrain = 0.0;
		    if (abs (tau0) > broken_limit) rVariables.EquivalentStrain = 0.0;
		    if (abs (tau1) > broken_limit) rVariables.EquivalentStrain = 0.0;
		}
	}
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void JointCohesionDriven3DLaw::CheckLoadingFunction(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    rVariables.LoadingFlag = false;
    rVariables.LoadingFunction = 0.0;

    if(rVariables.EquivalentStrain < mStateVariable)
    {
        rVariables.LoadingFlag = true;
        rVariables.LoadingFunction = 1.0;
    }
}

//----------------------------------------------------------------------------------------

void JointCohesionDriven3DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        // Unloading -> Tensile constitutive matrix
        if (mStateVariable == 1.0) // Unbroken joint
        {
			rConstitutiveMatrix(0,0) = rVariables.YieldStress;
			rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(2,2) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
			rConstitutiveMatrix(2,0) = 0.0;
			rConstitutiveMatrix(2,1) = 0.0;
			rConstitutiveMatrix(0,2) = 0.0;
			rConstitutiveMatrix(1,2) = 0.0;
		}

		if (mStateVariable == 0.0) // Broken joint
		{
            double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			rConstitutiveMatrix(0,0) = broken_YieldStress;
			rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(2,2) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
			rConstitutiveMatrix(2,0) = 0.0;
			rConstitutiveMatrix(2,1) = 0.0;
			rConstitutiveMatrix(0,2) = 0.0;
			rConstitutiveMatrix(1,2) = 0.0;
		}
    }

    else // Contact between interfaces
    {
        // Unloading -> Compresive constitutive matrix
        if (mStateVariable == 1.0) // Unbroken joint
        {
			rConstitutiveMatrix(0,0) = rVariables.YieldStress;
			rConstitutiveMatrix(1,1) = rVariables.YieldStress;
			rConstitutiveMatrix(2,2) = rVariables.YoungModulus;
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
			rConstitutiveMatrix(2,0) = 0.0;
			rConstitutiveMatrix(2,1) = 0.0;
			rConstitutiveMatrix(0,2) = 0.0;
			rConstitutiveMatrix(1,2) = 0.0;
		}

		if (mStateVariable == 0.0) // Broken joint
		{
			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			double shear_modulus_stress0 = fabs (StrainVector[0] / (2.0 * (1.0 + rVariables.PoissonCoefficient)));
			double shear_modulus_stress1 = fabs (StrainVector[1] / (2.0 * (1.0 + rVariables.PoissonCoefficient)));
			double friction_modulus_stress = fabs(rVariables.FrictionCoefficient * StrainVector[2]);

			if (shear_modulus_stress0 > friction_modulus_stress && shear_modulus_stress1 > friction_modulus_stress)
			{
				rConstitutiveMatrix(0,0) = broken_YieldStress;
				rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
				rConstitutiveMatrix(2,2) = rVariables.YoungModulus;
				rConstitutiveMatrix(0,1) = 0.0;
				rConstitutiveMatrix(1,0) = 0.0;
				rConstitutiveMatrix(2,0) = 0.0;
				rConstitutiveMatrix(2,1) = 0.0;

				const double eps = std::numeric_limits<double>::epsilon();

				if(StrainVector[0] > eps)
				{
					rConstitutiveMatrix(0,2) =-rVariables.YoungModulus * rVariables.FrictionCoefficient;
				}
				else if(StrainVector[0] < -eps)
				{
					rConstitutiveMatrix(0,2) = rVariables.YoungModulus * rVariables.FrictionCoefficient;
				}
				else
				{
					rConstitutiveMatrix(0,2) = 0.0;
				}

				if(StrainVector[1] > eps)
				{
					rConstitutiveMatrix(1,2) =-rVariables.YoungModulus * rVariables.FrictionCoefficient;
				}
				else if(StrainVector[1] < -eps)
				{
					rConstitutiveMatrix(1,2) = rVariables.YoungModulus * rVariables.FrictionCoefficient;
				}
				else
				{
					rConstitutiveMatrix(1,2) = 0.0;
				}
			}

			else
			{
				double shear_modulus = rVariables.YieldStress / (2.0 * (1.0 + rVariables.PoissonCoefficient));
				rConstitutiveMatrix(0,0) = broken_YieldStress + shear_modulus;
				rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
				rConstitutiveMatrix(2,2) = rVariables.YoungModulus;
				rConstitutiveMatrix(0,1) = 0.0;
				rConstitutiveMatrix(1,0) = 0.0;
				rConstitutiveMatrix(2,0) = 0.0;
				rConstitutiveMatrix(2,1) = 0.0;
				rConstitutiveMatrix(0,2) = 0.0;
				rConstitutiveMatrix(1,2) = 0.0;
			}
		}

    }
}

//----------------------------------------------------------------------------------------

void JointCohesionDriven3DLaw::ComputeStressVector(Vector& rStressVector,
                                                 ConstitutiveLawVariables& rVariables,
                                                 Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
		// Tensile stress
		if (mStateVariable == 1.0) // Unbroken joint
        {
			rStressVector[0] = rVariables.YieldStress * StrainVector[0];
			rStressVector[1] = rVariables.YieldStress * StrainVector[1];
			rStressVector[2] = rVariables.YieldStress * StrainVector[2];
		}

		if (mStateVariable == 0.0) // Broken joint
		{
			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			rStressVector[0] = broken_YieldStress * StrainVector[0];
			rStressVector[1] = broken_YieldStress * StrainVector[1];
			rStressVector[2] = broken_YieldStress * StrainVector[2];
		}


    }
    else // Contact between interfaces
    {
        // Note: StrainVector[1] < 0.0, rStressVector[1] < 0.0 -> Compresive stress
        if (mStateVariable == 1.0) // Unbroken joint
        {
			rStressVector[0] = rVariables.YieldStress * StrainVector[0];
			rStressVector[1] = rVariables.YieldStress * StrainVector[1];
			rStressVector[2] = rVariables.YoungModulus * StrainVector[2];
		}


        if (mStateVariable==0.0) // Broken joint
		{
			rStressVector[2] = rVariables.YoungModulus * StrainVector[2];

			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			double tangential_strain_vector_modulus = sqrt(StrainVector[0] * StrainVector[0] + StrainVector[1] * StrainVector[1]);
			const double shear_modulus =  rVariables.YieldStress / (2.0 * (1.0 + rVariables.PoissonCoefficient));

			double friction_stress0 = fabs(shear_modulus * StrainVector[0]);
			double friction_stress1 = fabs(shear_modulus * StrainVector[1]);
			double max_friction_stress = fabs(rVariables.FrictionCoefficient * rStressVector[2]);

			if (friction_stress0 > max_friction_stress) friction_stress0 = max_friction_stress;
			if (friction_stress1 > max_friction_stress) friction_stress1 = max_friction_stress;

			double friction_stress0_proj = fabs(friction_stress0 * StrainVector[0] / tangential_strain_vector_modulus);
			double friction_stress1_proj = fabs(friction_stress1 * StrainVector[1] / tangential_strain_vector_modulus);

			const double eps = std::numeric_limits<double>::epsilon();

			if(StrainVector[0] > eps)
			{
				rStressVector[0] = broken_YieldStress * StrainVector[0] + friction_stress0_proj;
			}
			else if(StrainVector[0] < -eps)
			{
				rStressVector[0] = broken_YieldStress * StrainVector[0] - friction_stress0_proj;
			}
			else
			{
				rStressVector[0] = 0.0;
			}

			if(StrainVector[1] > eps)
			{
				rStressVector[1] = broken_YieldStress * StrainVector[1] + friction_stress1_proj;
			}
			else if(StrainVector[1] < -eps)
			{
				rStressVector[1] = broken_YieldStress * StrainVector[1] - friction_stress1_proj;
			}
			else
			{
				rStressVector[1] = 0.0;
			}
		}
    }
}

} // Namespace Kratos
