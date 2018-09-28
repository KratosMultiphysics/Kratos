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
#include "custom_constitutive/hyperelastic_U_P_plane_strain_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticUPPlaneStrain2DLaw::HyperElasticUPPlaneStrain2DLaw()
    : HyperElasticUP3DLaw()
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticUPPlaneStrain2DLaw::HyperElasticUPPlaneStrain2DLaw(const HyperElasticUPPlaneStrain2DLaw& rOther)
    : HyperElasticUP3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticUPPlaneStrain2DLaw::Clone() const
{
    return Kratos::make_shared<HyperElasticUPPlaneStrain2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticUPPlaneStrain2DLaw::~HyperElasticUPPlaneStrain2DLaw()
{
}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticUPPlaneStrain2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );

}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticUPPlaneStrain2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
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


void HyperElasticUPPlaneStrain2DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
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


void HyperElasticUPPlaneStrain2DLaw::CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
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

void HyperElasticUPPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
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
