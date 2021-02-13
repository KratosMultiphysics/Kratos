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
#include "custom_constitutive/joint_stress_driven_3D_law.hpp"

namespace Kratos
{

int JointStressDriven3DLaw::Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    if(rMaterialProperties.Has(YOUNG_MODULUS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YOUNG_MODULUS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(POISSON_RATIO)) {
        KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] < -1.0) << "POISSON_RATIO has an invalid value lower than -1.0" << std::endl;
        KRATOS_ERROR_IF(rMaterialProperties[POISSON_RATIO] >= 0.5) << "POISSON_RATIO has an invalid value greater or equal to 0.5 " << std::endl;
    } else {
        KRATOS_ERROR << "POISSON_RATIO not defined" << std::endl;
    }

    if(rMaterialProperties.Has(MAX_COMPRESSIVE_STRESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[MAX_COMPRESSIVE_STRESS] < 0.0) << "MAX_COMPRESSIVE_STRESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "MAX_COMPRESSIVE_STRESS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(MAX_TENSILE_STRESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[MAX_TENSILE_STRESS] < 0.0) << "MAX_TENSILE_STRESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "MAX_TENSILE_STRESS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(FRICTION_COEFFICIENT)) {
        KRATOS_ERROR_IF(rMaterialProperties[FRICTION_COEFFICIENT] < 0.0) << "FRICTION_COEFFICIENT has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "FRICTION_COEFFICIENT not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------

void JointStressDriven3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable = 1.0;
}

//----------------------------------------------------------------------------------------

void JointStressDriven3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
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

void JointStressDriven3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)

{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    rVariables.YoungModulus = MaterialProperties[YOUNG_MODULUS];
    rVariables.PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    rVariables.MaxCompresiveStress = MaterialProperties[MAX_COMPRESSIVE_STRESS];
    rVariables.MaxTensileStress = MaterialProperties[MAX_TENSILE_STRESS];
    rVariables.YieldStress = rVariables.YoungModulus;
    rVariables.FrictionCoefficient = MaterialProperties[FRICTION_COEFFICIENT];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void JointStressDriven3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                     Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        rVariables.EquivalentStrain = 1.0;
        if (mStateVariable == 1.0)
        {
            if((sqrt(StrainVector[0] * StrainVector[0] + StrainVector[1] * StrainVector[1]) * rVariables.YieldStress) > rVariables.MaxTensileStress)
            {
			    rVariables.EquivalentStrain = 0.0;
			}
            if(fabs(rVariables.YieldStress * StrainVector[2]) > rVariables.MaxTensileStress)
            {
			    rVariables.EquivalentStrain = 0.0;
			}
        }
    }
    else // Contact between interfaces
    {
        rVariables.EquivalentStrain = 1.0;
        if (mStateVariable == 1.0)
        {
		    if((sqrt(StrainVector[0] * StrainVector[0] + StrainVector[1] * StrainVector[1]) * rVariables.YieldStress) > rVariables.MaxTensileStress)
            {
			    rVariables.EquivalentStrain = 0.0;
			}
            if(fabs(rVariables.YoungModulus * StrainVector[2]) > rVariables.MaxCompresiveStress)
            {
			    rVariables.EquivalentStrain = 0.0;
			}
		}
	}
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void JointStressDriven3DLaw::CheckLoadingFunction(ConstitutiveLawVariables& rVariables,
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

void JointStressDriven3DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
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

    }
}

//----------------------------------------------------------------------------------------

void JointStressDriven3DLaw::ComputeStressVector(Vector& rStressVector,
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
			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			rStressVector[2] = rVariables.YoungModulus * StrainVector[2];

			if (rVariables.FrictionCoefficient == 0.0)
			{
				rStressVector[0] = broken_YieldStress * StrainVector[0];
				rStressVector[1] = broken_YieldStress * StrainVector[1];
			}

			else {
				double tangential_strain_vector_modulus = sqrt(StrainVector[0] * StrainVector[0] + StrainVector[1] * StrainVector[1]);
				const double shear_modulus =  rVariables.YieldStress / (2.0 * (1.0 + rVariables.PoissonCoefficient));

				double friction_stress = fabs(shear_modulus * tangential_strain_vector_modulus);
				double max_friction_stress = fabs(rVariables.FrictionCoefficient * rStressVector[2]);

				if (friction_stress > max_friction_stress) friction_stress = max_friction_stress;

				double friction_stress0 = fabs(friction_stress * StrainVector[0] / tangential_strain_vector_modulus);
				double friction_stress1 = fabs(friction_stress * StrainVector[1] / tangential_strain_vector_modulus);

				const double eps = std::numeric_limits<double>::epsilon();

				if(StrainVector[0] > eps)
				{
					rStressVector[0] = broken_YieldStress * StrainVector[0] + friction_stress0;
				}
				else if(StrainVector[0] < -eps)
				{
					rStressVector[0] = broken_YieldStress * StrainVector[0] - friction_stress0;
				}
				else
				{
					rStressVector[0] = 0.0;
				}

				if(StrainVector[1] > eps)
				{
					rStressVector[1] = broken_YieldStress * StrainVector[1] + friction_stress1;
				}
				else if(StrainVector[1] < -eps)
				{
					rStressVector[1] = broken_YieldStress * StrainVector[1] - friction_stress1;
				}
				else
				{
					rStressVector[1] = 0.0;
				}
			}
		}
    }
}

} // Namespace Kratos
