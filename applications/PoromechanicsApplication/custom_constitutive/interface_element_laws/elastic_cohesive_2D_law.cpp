//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana and Danilo Cavalcanti
//

// Application includes
#include "custom_constitutive/interface_element_laws/elastic_cohesive_2D_law.hpp"

namespace Kratos
{

void ElasticCohesive2DLaw::GetLawFeatures(Features& rFeatures)
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

//----------------------------------------------------------------------------------------

void ElasticCohesive2DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
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
    noalias(rConstitutiveMatrix) = ZeroMatrix(2,2);
    rConstitutiveMatrix(0,0) = rVariables.ShearStiffness;
    rConstitutiveMatrix(1,1) = rVariables.NormalStiffness * cp;

}

//----------------------------------------------------------------------------------------

void ElasticCohesive2DLaw::ComputeStressVector(Vector& rStressVector,
                                                ConstitutiveLawVariables& rVariables,
                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    double cp = 1.0;

    // Penalization coefficient, in case it is a compression
    if(StrainVector[1] < 0.0){
        cp = rVariables.PenaltyStiffness;
    }

    rStressVector[0] = StrainVector[0] * rVariables.ShearStiffness;
    rStressVector[1] = StrainVector[1] * rVariables.NormalStiffness * cp;
    
}

} // Namespace Kratos
