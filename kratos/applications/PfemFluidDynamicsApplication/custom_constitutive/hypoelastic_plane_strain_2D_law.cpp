//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Last modified by:    $Author:                   AFranci $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"

#include "solid_mechanics_application.h"

#include "pfem_fluid_dynamics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HypoElasticPlaneStrain2DLaw::HypoElasticPlaneStrain2DLaw()
    : HypoElastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HypoElasticPlaneStrain2DLaw::HypoElasticPlaneStrain2DLaw(const HypoElasticPlaneStrain2DLaw& rOther)
    : HypoElastic3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HypoElasticPlaneStrain2DLaw::Clone() const
{
    HypoElasticPlaneStrain2DLaw::Pointer p_clone(new HypoElasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HypoElasticPlaneStrain2DLaw::~HypoElasticPlaneStrain2DLaw()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HypoElasticPlaneStrain2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );

}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HypoElasticPlaneStrain2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen ( 2 , 2 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = -InverseLeftCauchyGreen( 0, 1 );


}



//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************

void HypoElasticPlaneStrain2DLaw::CalculateConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
        Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();


    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables,
                                          this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
        }

    }


}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HypoElasticPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
	
	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}


} // Namespace Kratos
