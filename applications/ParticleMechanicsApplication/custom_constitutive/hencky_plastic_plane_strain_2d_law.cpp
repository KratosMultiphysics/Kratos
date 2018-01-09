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

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_plane_strain_2d_law.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrain2DLaw::HenckyElasticPlasticPlaneStrain2DLaw()
    : HenckyElasticPlastic3DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrain2DLaw::HenckyElasticPlasticPlaneStrain2DLaw(FlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HenckyElasticPlastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrain2DLaw::HenckyElasticPlasticPlaneStrain2DLaw(const HenckyElasticPlasticPlaneStrain2DLaw& rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlasticPlaneStrain2DLaw::Clone() const
{
    HenckyElasticPlasticPlaneStrain2DLaw::Pointer p_clone(new HenckyElasticPlasticPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticPlaneStrain2DLaw::~HenckyElasticPlasticPlaneStrain2DLaw()
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticPlaneStrain2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );

}


//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticPlaneStrain2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen ( rLeftCauchyGreen.size1(), rLeftCauchyGreen.size2() );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = -InverseLeftCauchyGreen( 0, 1 );

    //std::cout<<" in CalculateAlmansiStrain "<<std::endl;


}


Vector HenckyElasticPlasticPlaneStrain2DLaw::SetStressMatrixToAppropiateVectorDimension(Vector& rStressVector, const Matrix& rStressMatrix)
{
    //std::cout<<" rStressMatrix " <<rStressMatrix<<std::endl;
    rStressVector = MathUtils<double>::StressTensorToVector( rStressMatrix, rStressVector.size() );
    //if(rStressVector.size() == 6)
    //{
    //rStressVector = ZeroVector(6);
    //rStressVector(0) = rStressMatrix(0,0);
    //rStressVector(1) = rStressMatrix(1,1);
    //rStressVector(2) = rStressMatrix(2,2);
    //rStressVector(3) = rStressMatrix(0,1);
    //}
    //else
    //{
    //rStressVector = ZeroVector(3);
    //rStressVector(0) = rStressMatrix(0,0);
    //rStressVector(1) = rStressMatrix(1,1);
    //rStressVector(2) = rStressMatrix(0,1);
    //}
    //std::cout<<" rStressVector " <<rStressVector<<std::endl;
    return rStressVector;

}


Matrix HenckyElasticPlasticPlaneStrain2DLaw::SetConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix, const Matrix& rElastoPlasticTangentMatrix)
{
    if(rConstitutiveMatrix.size1() == 6)
    {
        rConstitutiveMatrix = ZeroMatrix(6,6);
        rConstitutiveMatrix = rElastoPlasticTangentMatrix;
        //std::cout<<" rConstitutiveMatrix " <<rConstitutiveMatrix<<std::endl;
    }
    else
    {
        rConstitutiveMatrix = ZeroMatrix(3,3);
        //ilaria: it has been modified for the use of II mixed formulation

        //std::cout<<" rElastoPlasticTangentMatrix " <<rElastoPlasticTangentMatrix<<std::endl;
        rConstitutiveMatrix(0, 0) = rElastoPlasticTangentMatrix(0, 0);
        rConstitutiveMatrix(0, 1) = rElastoPlasticTangentMatrix(0, 1);
        rConstitutiveMatrix(1, 0) = rElastoPlasticTangentMatrix(1, 0);
        rConstitutiveMatrix(1, 1) = rElastoPlasticTangentMatrix(1, 1);


        rConstitutiveMatrix(2, 0) = rElastoPlasticTangentMatrix(3, 0);
        rConstitutiveMatrix(2, 1) = rElastoPlasticTangentMatrix(3, 1);
        rConstitutiveMatrix(2, 2) = rElastoPlasticTangentMatrix(3, 3);


        rConstitutiveMatrix(0, 2) = rElastoPlasticTangentMatrix(0, 3);
        rConstitutiveMatrix(1, 2) = rElastoPlasticTangentMatrix(1, 3);
        //std::cout<<" rConstitutiveMatrix " <<rConstitutiveMatrix<<std::endl;
    }
    return rConstitutiveMatrix;

}



void HenckyElasticPlasticPlaneStrain2DLaw::CalculateHenckyMainStrain(const Matrix& rCauchyGreenMatrix,
        MPMFlowRule::RadialReturnVariables& rReturnMappingVariables,
        Vector& rMainStrain)
{
    Matrix Auxiliar = ZeroMatrix(3,3);
    Auxiliar(0,0) = rCauchyGreenMatrix(0,0);
    Auxiliar(1,1) = rCauchyGreenMatrix(1,1);
    Auxiliar(0,1) = rCauchyGreenMatrix(0,1);
    Auxiliar(1,0) = rCauchyGreenMatrix(1,0);
    Auxiliar(2,2) = 1.0;
    Matrix AuxEigenVectors = ZeroMatrix(3,3);
    Vector AuxEigenValues  = ZeroVector(3);
    SolidMechanicsMathUtilities<double>::EigenVectors(Auxiliar, AuxEigenVectors, AuxEigenValues);


    Matrix EigenVectors = ZeroMatrix(3,3);
    EigenVectors(0,0) = AuxEigenVectors(0,0);
    EigenVectors(1,0) = AuxEigenVectors(1,0);
    EigenVectors(1,1) = AuxEigenVectors(1,1);
    EigenVectors(0,1) = AuxEigenVectors(0,1);
    // Positions known to be zero
    EigenVectors(0,2) = 0.0;
    EigenVectors(1,2) = 0.0;
    EigenVectors(2,0) = 0.0;
    EigenVectors(2,1) = 0.0;
    EigenVectors(2,2) = 1.0;

    rReturnMappingVariables.MainDirections     = EigenVectors;


    Vector TrialEigenValues = ZeroVector(3);
    TrialEigenValues(0) = AuxEigenValues(0);
    TrialEigenValues(1) = AuxEigenValues(1);
    TrialEigenValues(2) =  rCauchyGreenMatrix(2,2);

    for (unsigned int i = 0; i<3; i++)
        rMainStrain(i) = 0.50*std::log(TrialEigenValues(i));

    //std::cout<<" rCauchyGreenMatrix "<<rCauchyGreenMatrix<<std::endl;
    //std::cout<<" EigenVectors "<<EigenVectors<<std::endl;
    //std::cout<<" AuxEigenValues "<<AuxEigenValues<<std::endl;
    //std::cout<<" rMainStrain "<<rMainStrain<<std::endl;
}







} // Namespace Kratos

