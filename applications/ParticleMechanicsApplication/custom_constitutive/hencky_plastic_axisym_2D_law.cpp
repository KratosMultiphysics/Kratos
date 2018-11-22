//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_axisym_2d_law.hpp"
#include "custom_utilities/particle_mechanics_math_utilities.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::HenckyElasticPlasticAxisym2DLaw()
    : HenckyElasticPlastic3DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::HenckyElasticPlasticAxisym2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HenckyElasticPlastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::HenckyElasticPlasticAxisym2DLaw(const HenckyElasticPlasticAxisym2DLaw& rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlasticAxisym2DLaw::Clone() const
{
    HenckyElasticPlasticAxisym2DLaw::Pointer p_clone(new HenckyElasticPlasticAxisym2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticAxisym2DLaw::~HenckyElasticPlasticAxisym2DLaw()
{
}

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticAxisym2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    // E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 );

}


//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticAxisym2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
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


Matrix HenckyElasticPlasticAxisym2DLaw::SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix, const Matrix& rElastoPlasticTangentMatrix)
{
    if(rConstitutiveMatrix.size1() == 4)
    {
        rConstitutiveMatrix = ZeroMatrix(4,4);
        for(unsigned int i = 0; i<4; i++){
            for(unsigned int j = 0; j<4; j++){
                rConstitutiveMatrix(i,j) = rElastoPlasticTangentMatrix(i,j);
            }
	    }
    }
    else if(rConstitutiveMatrix.size1() == 6)
    {
        rConstitutiveMatrix = ZeroMatrix(6,6);
        rConstitutiveMatrix = rElastoPlasticTangentMatrix;
    }
    else if(rConstitutiveMatrix.size1() == 3)
    {
        //NOTE: It has been modified for the use of II mixed formulation
        rConstitutiveMatrix = ZeroMatrix(3,3);

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

} // Namespace Kratos

