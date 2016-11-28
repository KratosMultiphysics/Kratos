//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_U_P_axisym_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticUPAxisym2DLaw::HyperElasticUPAxisym2DLaw()
    : HyperElasticUP3DLaw()
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticUPAxisym2DLaw::HyperElasticUPAxisym2DLaw(const HyperElasticUPAxisym2DLaw& rOther)
    : HyperElasticUP3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticUPAxisym2DLaw::Clone() const
{
    HyperElasticUPAxisym2DLaw::Pointer p_clone(new HyperElasticUPAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticUPAxisym2DLaw::~HyperElasticUPAxisym2DLaw()
{
}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticUPAxisym2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 );

}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticUPAxisym2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen ( 3 , 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 2, 2 ) );
    rStrainVector[3] = -InverseLeftCauchyGreen( 0, 1 );


}


//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************


void HyperElasticUPAxisym2DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
        const Matrix & rIsoStressMatrix,
        Matrix& rConstitutiveMatrix)
{
    rConstitutiveMatrix.clear();

    for(unsigned int i=0; i<4; i++)
    {
        for(unsigned int j=0; j<4; j++)
        {
            rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rIsoStressMatrix,
                                          this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
        }

    }


}

//***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
//************************************************************************************


void HyperElasticUPAxisym2DLaw::CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
								       Matrix& rConstitutiveMatrix)

{
    rConstitutiveMatrix.clear();

    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors = this->CalculateVolumetricPressureFactors( rElasticVariables, Factors );

    for(unsigned int i=0; i<4; i++)
    {
        for(unsigned int j=0; j<4; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, Factors,
                                          this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
        }

    }


}



//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticUPAxisym2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
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
