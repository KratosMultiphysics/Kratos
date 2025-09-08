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
//  References:      This class is adapted from applications/MPMApplication/custom_constitutive/hyperelastic_plane_strain_2D_law.cpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/displacement_newtonian_fluid_plane_strain_2D_law.hpp"

#include "mpm_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluidPlaneStrain2DLaw::DispNewtonianFluidPlaneStrain2DLaw()
    : DispNewtonianFluid3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

DispNewtonianFluidPlaneStrain2DLaw::DispNewtonianFluidPlaneStrain2DLaw(const DispNewtonianFluidPlaneStrain2DLaw& rOther)
    : DispNewtonianFluid3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer DispNewtonianFluidPlaneStrain2DLaw::Clone() const
{
    return Kratos::make_shared<DispNewtonianFluidPlaneStrain2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluidPlaneStrain2DLaw::~DispNewtonianFluidPlaneStrain2DLaw()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************COMPUTE TOTAL STRAIN VECTOR**********************************
//************************************************************************************

void DispNewtonianFluidPlaneStrain2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
						Vector& rStrainVector )
{
    // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))
    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen = ZeroMatrix( 2, 2 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy

}


//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************

void DispNewtonianFluidPlaneStrain2DLaw::CalculateConstitutiveMatrix (const MaterialResponseVariables& rViscousVariables, Matrix& rConstitutiveMatrix)
{
    rConstitutiveMatrix.clear();

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ), rViscousVariables,
                                          this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1], this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
        }
    }

}



//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void DispNewtonianFluidPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
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
