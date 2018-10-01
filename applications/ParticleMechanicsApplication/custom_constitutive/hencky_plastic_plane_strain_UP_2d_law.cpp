//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//


// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_plane_strain_UP_2d_law.hpp"
//#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrainUP2DLaw::HenckyElasticPlasticPlaneStrainUP2DLaw()
    : HenckyElasticPlasticUP3DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrainUP2DLaw::HenckyElasticPlasticPlaneStrainUP2DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HenckyElasticPlasticUP3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrainUP2DLaw::HenckyElasticPlasticPlaneStrainUP2DLaw(const HenckyElasticPlasticPlaneStrainUP2DLaw& rOther)
    : HenckyElasticPlasticUP3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlasticPlaneStrainUP2DLaw::Clone() const
{
    HenckyElasticPlasticPlaneStrainUP2DLaw::Pointer p_clone(new HenckyElasticPlasticPlaneStrainUP2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrainUP2DLaw::~HenckyElasticPlasticPlaneStrainUP2DLaw()
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticPlaneStrainUP2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{
    // E = 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 ); // xy
}


//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticPlaneStrainUP2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // E = 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen ( rLeftCauchyGreen.size1(), rLeftCauchyGreen.size2() );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = -InverseLeftCauchyGreen( 0, 1 );

}


Vector HenckyElasticPlasticPlaneStrainUP2DLaw::SetStressMatrixToAppropiateVectorDimension(Vector& rStressVector, const Matrix& rStressMatrix)
{
    rStressVector = MathUtils<double>::StressTensorToVector( rStressMatrix, rStressVector.size() );
    return rStressVector;
}


Matrix HenckyElasticPlasticPlaneStrainUP2DLaw::SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix, const Matrix& rElastoPlasticTangentMatrix)
{
    // It has been modified for the use of mixed formulation
    if(rConstitutiveMatrix.size1() == 6)
    {
        rConstitutiveMatrix = rElastoPlasticTangentMatrix;
    }
    else
    {
        rConstitutiveMatrix(0, 0) = rElastoPlasticTangentMatrix(0, 0);
        rConstitutiveMatrix(0, 1) = rElastoPlasticTangentMatrix(0, 1);
        rConstitutiveMatrix(1, 0) = rElastoPlasticTangentMatrix(1, 0);
        rConstitutiveMatrix(1, 1) = rElastoPlasticTangentMatrix(1, 1);

        rConstitutiveMatrix(2, 0) = rElastoPlasticTangentMatrix(3, 0);
        rConstitutiveMatrix(2, 1) = rElastoPlasticTangentMatrix(3, 1);
        rConstitutiveMatrix(2, 2) = rElastoPlasticTangentMatrix(3, 3);

        rConstitutiveMatrix(0, 2) = rElastoPlasticTangentMatrix(0, 3);
        rConstitutiveMatrix(1, 2) = rElastoPlasticTangentMatrix(1, 3);
    }
    return rConstitutiveMatrix;

}

void HenckyElasticPlasticPlaneStrainUP2DLaw::GetLawFeatures(Features& rFeatures)
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

