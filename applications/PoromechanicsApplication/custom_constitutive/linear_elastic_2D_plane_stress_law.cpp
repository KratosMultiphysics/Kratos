//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Project includes
#include "custom_constitutive/linear_elastic_2D_plane_stress_law.hpp"

namespace Kratos
{

//Default Constructor
LinearElastic2DPlaneStressLaw::LinearElastic2DPlaneStressLaw() : LinearElastic3DLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
LinearElastic2DPlaneStressLaw::LinearElastic2DPlaneStressLaw(const LinearElastic2DPlaneStressLaw& rOther) : LinearElastic3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
LinearElastic2DPlaneStressLaw::~LinearElastic2DPlaneStressLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic2DPlaneStressLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRESS_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 2;
    
	//Set the strain size
	rFeatures.mStrainSize = 3;
}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer LinearElastic2DPlaneStressLaw::Clone() const
{
    LinearElastic2DPlaneStressLaw::Pointer p_clone(new LinearElastic2DPlaneStressLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void LinearElastic2DPlaneStressLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,const double& rYoungModulus,const double& rPoissonCoefficient )
{
    rConstitutiveMatrix.clear();

    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus)/(1.0-rPoissonCoefficient*rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 )*(1-rPoissonCoefficient)*0.5;

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient;
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );
}

} // Namespace Kratos
