//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/hyperelastic_U_P_plane_strain_2D_law.hpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_plane_strain_UP_2D_law.hpp"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlaneStrainUP2DLaw::HyperElasticPlaneStrainUP2DLaw()
    : HyperElasticUP3DLaw()
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlaneStrainUP2DLaw::HyperElasticPlaneStrainUP2DLaw(const HyperElasticPlaneStrainUP2DLaw& rOther)
    : HyperElasticUP3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlaneStrainUP2DLaw::Clone() const
{
    return Kratos::make_shared<HyperElasticPlaneStrainUP2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlaneStrainUP2DLaw::~HyperElasticPlaneStrainUP2DLaw()
{
}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticPlaneStrainUP2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );

}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticPlaneStrainUP2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
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


//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************


void HyperElasticPlaneStrainUP2DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
									   const Matrix & rIsoStressMatrix,
									   Matrix& rConstitutiveMatrix)
{
    rConstitutiveMatrix.clear();

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rIsoStressMatrix,
                                          this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
        }

    }


}

//***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
//************************************************************************************


void HyperElasticPlaneStrainUP2DLaw::CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
									    Matrix& rConstitutiveMatrix)

{
    rConstitutiveMatrix.clear();

    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors = this->CalculateVolumetricPressureFactors( rElasticVariables, Factors );

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, Factors,
                                          this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
        }

    }


}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticPlaneStrainUP2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );
	rFeatures.mOptions.Set( U_P_LAW );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}


} // Namespace Kratos
