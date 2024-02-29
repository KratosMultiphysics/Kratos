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
#include "custom_constitutive/interface_element_laws/elastoplastic_mod_mohr_coulomb_cohesive_2D_law.hpp"

namespace Kratos
{

void ElastoPlasticModMohrCoulombCohesive2DLaw::GetLawFeatures(Features& rFeatures)
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

void ElastoPlasticModMohrCoulombCohesive2DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mPlasticStrainVector.resize(2);
    mOldPlasticStrainVector.resize(2);
    mPlasticStrainVector[0]    = 0.0;
    mPlasticStrainVector[1]    = 0.0;
    mOldPlasticStrainVector[0] = 0.0;
    mOldPlasticStrainVector[1] = 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method returns the resultant shear component of the stress vector 

double ElastoPlasticModMohrCoulombCohesive2DLaw::GetShearResultantStressVector(Vector& StressVector)
{
    // Compute the shear resultant of the stress vector
    return StressVector[0];
}

//----------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive2DLaw::GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix,
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
