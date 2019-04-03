//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlaneStress2DLaw::LinearElasticPlaneStress2DLaw()
    : LinearElasticPlaneStrain2DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearElasticPlaneStress2DLaw::LinearElasticPlaneStress2DLaw(const LinearElasticPlaneStress2DLaw& rOther)
    : LinearElasticPlaneStrain2DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearElasticPlaneStress2DLaw::Clone() const
{
    return Kratos::make_shared<LinearElasticPlaneStress2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlaneStress2DLaw::~LinearElasticPlaneStress2DLaw()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************



//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************


void LinearElasticPlaneStress2DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
        const double &rYoungModulus,
        const double &rPoissonCoefficient )
{
    rConstitutiveMatrix.clear();

    // Plane stress constitutive matrix
    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus)/(1.0-rPoissonCoefficient*rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-rPoissonCoefficient)*0.5;

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient;
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearElasticPlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRESS_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}


} // Namespace Kratos
