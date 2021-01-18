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
#include "custom_constitutive/joint_stress_driven_2D_law.hpp"

namespace Kratos
{

void JointStressDriven2DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                     Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        rVariables.EquivalentStrain = 1.0;
        if (mStateVariable == 1.0)
        {
            if(fabs(StrainVector[0] * rVariables.YieldStress) > rVariables.MaxTensileStress)
            {
			    rVariables.EquivalentStrain = 0.0;
			}
            if(fabs(rVariables.YieldStress * StrainVector[1]) > rVariables.MaxTensileStress)
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
		    if(fabs(StrainVector[0] * rVariables.YieldStress) > rVariables.MaxTensileStress)
            {
			    rVariables.EquivalentStrain = 0.0;
			}
            if(fabs(rVariables.YoungModulus * StrainVector[1]) > rVariables.MaxCompresiveStress)
            {
			    rVariables.EquivalentStrain = 0.0;
			}
		}
	}
}

//----------------------------------------------------------------------------------------

void JointStressDriven2DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                       ConstitutiveLawVariables& rVariables,
                                                       Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) ) // No contact between interfaces
    {
        // Tensile constitutive matrix
        if (mStateVariable == 1.0) // Unbroken joint
        {
			rConstitutiveMatrix(0,0) = rVariables.YieldStress;
			rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
		}

		if (mStateVariable == 0.0) // Broken joint
		{
			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			rConstitutiveMatrix(0,0) = broken_YieldStress;
			rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
			rConstitutiveMatrix(0,1) = 0.0;
			rConstitutiveMatrix(1,0) = 0.0;
		}
    }

    else // Contact between interfaces
    {
        // Compresive constitutive matrix
        if (mStateVariable == 1.0) // Unbroken joint
		{
			rConstitutiveMatrix(0,0) = rVariables.YieldStress;
		}

		if (mStateVariable == 0.0) // Broken joint
		{
			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			rConstitutiveMatrix(0,0) = broken_YieldStress;
		}

        rConstitutiveMatrix(1,1) = rVariables.YoungModulus;

        const double eps = std::numeric_limits<double>::epsilon();

        if(StrainVector[0] > eps)
        {
            rConstitutiveMatrix(0,1) =-rVariables.YoungModulus * rVariables.FrictionCoefficient;
        }
        else if(StrainVector[0] < -eps)
        {
            rConstitutiveMatrix(0,1) = rVariables.YoungModulus * rVariables.FrictionCoefficient;
        }
        else
        {
            rConstitutiveMatrix(0,1) = 0.0;
        }

        rConstitutiveMatrix(1,0) = 0.0;
    }
}

//----------------------------------------------------------------------------------------

void JointStressDriven2DLaw::ComputeStressVector(Vector& rStressVector,
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
		}

		if (mStateVariable == 0.0) // Broken joint
		{
			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			rStressVector[0] = broken_YieldStress * StrainVector[0];
			rStressVector[1] = broken_YieldStress * StrainVector[1];
		}
    }

    else // Contact between interfaces
    {
        // Note: StrainVector[1] < 0.0, rStressVector[1] < 0.0 -> Compresive stress

        if (mStateVariable==1.0) // Unbroken joint
		{
			rStressVector[0] = rVariables.YieldStress * StrainVector[0];
			rStressVector[1] = rVariables.YoungModulus * StrainVector[1];
		}

        if (mStateVariable==0.0) // Broken joint
		{
			double broken_YieldStress = rVariables.YoungModulus * 1.0e-9;
			rStressVector[1] = broken_YieldStress * StrainVector[1];

			if (rVariables.FrictionCoefficient == 0.0) rStressVector[0] = broken_YieldStress * StrainVector[0];

			else {
				const double shear_modulus = rVariables.YieldStress / (2.0 * (1.0 + rVariables.PoissonCoefficient));
				double friction_stress = fabs(shear_modulus * StrainVector[0]);
				double max_friction_stress = fabs(rVariables.FrictionCoefficient * rStressVector[1]);
				if (friction_stress > max_friction_stress) friction_stress = max_friction_stress;

				const double eps = std::numeric_limits<double>::epsilon();

				if(StrainVector[0] > eps)
				{
					rStressVector[0] = broken_YieldStress * StrainVector[0] + friction_stress;
				}
				else if(StrainVector[0] < -eps)
				{
					rStressVector[0] = broken_YieldStress * StrainVector[0] - friction_stress;
				}
				else
				{
					rStressVector[0] = 0.0;
				}
			}
		}
    }
}

} // Namespace Kratos
