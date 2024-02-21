//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Danilo Cavalcanti and Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/interface_element_laws/isotropic_damage_cohesive_2D_law.hpp"

namespace Kratos
{

void IsotropicDamageCohesive2DLaw::GetLawFeatures(Features& rFeatures)
{
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
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the current state variables, current equivalent strain, and the vector with the 
// derivatives of the equivalent strain wrt to the components of the strain vector

void IsotropicDamageCohesive2DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    // Maximum shear resultant and normal strains in the loading history
    mStateVariable[0] = std::max(mOldStateVariable[0],StrainVector[0]);
    mStateVariable[1] = std::max(mOldStateVariable[1],StrainVector[1]);

    // Compute the current equivalent strain
    rVariables.EquivalentStrain = mStateVariable[1] + rVariables.BetaEqStrainShearFactor * mStateVariable[0];

    // Compute the equivalent strain associated with the last converged step (used to evaluate the loading function)
    rVariables.OldEquivalentStrain = mOldStateVariable[1] + rVariables.BetaEqStrainShearFactor * mOldStateVariable[0];

    // Compute the vector with the derivatives of the equivalent strain wrt to the components of the strain vector
    const double sign = (StrainVector[0] < 0.0) ? -1.0 : 1.0;
    rVariables.DerivativeEquivalentStrain[0] = rVariables.BetaEqStrainShearFactor * sign;
    rVariables.DerivativeEquivalentStrain[1] = 1.0;
}

//----------------------------------------------------------------------------------------

void IsotropicDamageCohesive2DLaw::GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    double cp = 1.0;
    // Penalization coefficient, in case it is a compression
    if(StrainVector[1] < 0.0){
        cp = rVariables.PenaltyStiffness;
    }

    // Fill the constitutive matrix
    noalias(rElasticConstitutiveMatrix) = ZeroMatrix(2,2);
    rElasticConstitutiveMatrix(0,0) = rVariables.ShearStiffness;
    rElasticConstitutiveMatrix(1,1) = rVariables.NormalStiffness * cp;
}

} // Namespace Kratos
