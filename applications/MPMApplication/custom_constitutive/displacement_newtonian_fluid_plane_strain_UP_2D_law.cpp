//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Contri Alessandro
//
//  References:      This class is adapted from applications/ParticleMechanicsApplication/custom_constitutive/hyperelastic_UP_plane_strain_2D_law.hpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/displacement_newtonian_fluid_plane_strain_UP_2D_law.hpp"
#include "mpm_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluidPlaneStrainUP2DLaw::DispNewtonianFluidPlaneStrainUP2DLaw()
    : DispNewtonianFluidUP3DLaw()
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

DispNewtonianFluidPlaneStrainUP2DLaw::DispNewtonianFluidPlaneStrainUP2DLaw(const DispNewtonianFluidPlaneStrainUP2DLaw& rOther)
    : DispNewtonianFluidUP3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer DispNewtonianFluidPlaneStrainUP2DLaw::Clone() const
{
    return Kratos::make_shared<DispNewtonianFluidPlaneStrainUP2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluidPlaneStrainUP2DLaw::~DispNewtonianFluidPlaneStrainUP2DLaw()
{
}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void DispNewtonianFluidPlaneStrainUP2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen = ZeroMatrix( 2 , 2 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = -InverseLeftCauchyGreen( 0, 1 );


}


//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************


void DispNewtonianFluidPlaneStrainUP2DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rViscousVariables,
									   const Matrix & rIsoStressMatrix,
									   Matrix& rConstitutiveMatrix)
{
    rConstitutiveMatrix.clear();

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rViscousVariables, rIsoStressMatrix,
                                          this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
        }

    }


}

//***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
//************************************************************************************


void DispNewtonianFluidPlaneStrainUP2DLaw::CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables & rViscousVariables,
									    Matrix& rConstitutiveMatrix)

{
    rConstitutiveMatrix.clear();

    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors = this->CalculateVolumetricPressureFactors( rViscousVariables, Factors );

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rViscousVariables, Factors,
                                          this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
        }

    }


}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void DispNewtonianFluidPlaneStrainUP2DLaw::GetLawFeatures(Features& rFeatures)
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
