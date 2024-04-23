//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/dam_joint_2D_law.hpp"

namespace Kratos
{

void DamJoint2DLaw::GetLawFeatures(Features& rFeatures)
{
    KRATOS_TRY

    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	//rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 2;

	//Set the strain size
	rFeatures.mStrainSize = 2;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void DamJoint2DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                            Parameters& rValues)
{}

//----------------------------------------------------------------------------------------

void DamJoint2DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                              ConstitutiveLawVariables& rVariables,
                                              Parameters& rValues)
{
    KRATOS_TRY

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY)) // No contact between interfaces
    {
        rConstitutiveMatrix(0,0) = rVariables.YieldStress;
        rConstitutiveMatrix(1,1) = rVariables.YieldStress;

        rConstitutiveMatrix(0,1) = 0.0;
        rConstitutiveMatrix(1,0) = 0.0;
    }

    else // Contact between interfaces
    {
        rConstitutiveMatrix(0,0) = rVariables.YieldStress;
        rConstitutiveMatrix(1,1) = rVariables.YoungModulus;

        rConstitutiveMatrix(0,1) = 0.0;
        rConstitutiveMatrix(1,0) = 0.0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void DamJoint2DLaw::ComputeStressVector(Vector& rStressVector,
                                        ConstitutiveLawVariables& rVariables,
                                        Parameters& rValues)
{
    KRATOS_TRY

    const Vector& StrainVector = rValues.GetStrainVector();

    if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY)) // No contact between interfaces
    {
        rStressVector[0] = rVariables.YieldStress * StrainVector[0];
        rStressVector[1] = rVariables.YieldStress * StrainVector[1];
    }
    else
    {
        rStressVector[0] = rVariables.YieldStress * StrainVector[0];
        rStressVector[1] = rVariables.YoungModulus * StrainVector[1];
    }

    KRATOS_CATCH("")
}

} // Namespace Kratos
