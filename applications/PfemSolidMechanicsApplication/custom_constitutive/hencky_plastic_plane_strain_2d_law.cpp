//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_plane_strain_2d_law.hpp"

#include "solid_mechanics_application.h"

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

HenckyElasticPlasticPlaneStrain2DLaw::HenckyElasticPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
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
    Matrix InverseLeftCauchyGreen ( rLeftCauchyGreen.size1() , rLeftCauchyGreen.size2() );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = -InverseLeftCauchyGreen( 0, 1 );


}


//*********************** COMPUTE THE FIRST TERM OF THE CONSTITUTIVE MATRIX **********
//************************************************************************************

void HenckyElasticPlasticPlaneStrain2DLaw::CalculateEigenValuesConstitutiveMatrix(const Matrix & rPrincipalTangent, const Matrix & rEigenVectors, Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();


    static const unsigned int msIndexVoigt2D [3][2] = { {0, 0}, {1, 1}, {0, 1} };

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = EigenValuesConstitutiveComponent(rConstitutiveMatrix( i, j ), 
					  rPrincipalTangent, rEigenVectors, 
                                          msIndexVoigt2D[i][0], msIndexVoigt2D[i][1], msIndexVoigt2D[j][0], msIndexVoigt2D[j][1]);
        }

    }

}

//*********************** COMPUTE THE SECOND TERM OF THE CONSTITUTIVE MATRIX**********
//************************************************************************************
void HenckyElasticPlasticPlaneStrain2DLaw::CalculateEigenVectorsConstitutiveMatrix(const Matrix& rEigenVectors, const Vector& rEigenValues, const Matrix& rStressMatrix, Matrix& rConstitutiveMatrix)
{
    Vector PrincipalStress;
    PrincipalStress = this->GetStressVectorFromMatrix(rStressMatrix, PrincipalStress);
    
    rConstitutiveMatrix.clear();


    static const unsigned int msIndexVoigt2D [3][2] = { {0, 0}, {1, 1}, {0, 1} };

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = EigenVectorsConstitutiveComponent(rConstitutiveMatrix( i, j ), 
					  rEigenVectors, rEigenValues, PrincipalStress,
                                          msIndexVoigt2D[i][0], msIndexVoigt2D[i][1], msIndexVoigt2D[j][0], msIndexVoigt2D[j][1]);
        }

    }

}




void HenckyElasticPlasticPlaneStrain2DLaw::CalculatePrincipalAxisHenckyTrial(const Matrix& rCauchyGreenMatrix, FlowRule::RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStrain)
{
    Matrix Auxiliar = ZeroMatrix(3);
    Auxiliar(0,0) = rCauchyGreenMatrix(0,0);
    Auxiliar(1,1) = rCauchyGreenMatrix(1,1);
    Auxiliar(0,1) = rCauchyGreenMatrix(0,1);
    Auxiliar(1,0) = rCauchyGreenMatrix(1,0);
    Auxiliar(2,2) = 1.0; 
    Matrix AuxEigenVectors = ZeroMatrix(3);
    Vector AuxEigenValues  = ZeroVector(3);
    SD_MathUtils<double>::EigenVectors(Auxiliar, AuxEigenVectors, AuxEigenValues);

    
    rReturnMappingVariables.EigenVectors = ZeroMatrix(3);
    rReturnMappingVariables.EigenVectors(0,0) = AuxEigenVectors(0,0);
    rReturnMappingVariables.EigenVectors(1,0) = AuxEigenVectors(1,0);
    rReturnMappingVariables.EigenVectors(1,1) = AuxEigenVectors(1,1);
    rReturnMappingVariables.EigenVectors(0,1) = AuxEigenVectors(0,1);
    // Positions known to be zero
    rReturnMappingVariables.EigenVectors(0,2) = 0.0;
    rReturnMappingVariables.EigenVectors(1,2) = 0.0;
    rReturnMappingVariables.EigenVectors(2,0) = 0.0;
    rReturnMappingVariables.EigenVectors(2,1) = 0.0;
    rReturnMappingVariables.EigenVectors(2,2) = 1.0;

    rReturnMappingVariables.TrialEigenValues = ZeroVector(3);
    rReturnMappingVariables.TrialEigenValues(0) = AuxEigenValues(0);
    rReturnMappingVariables.TrialEigenValues(1) = AuxEigenValues(1);
    rReturnMappingVariables.TrialEigenValues(2) = 1.0; // rCauchyGreenMatrix(2,2);

    for (unsigned int i = 0; i<3; i++)
           rPrincipalStrain(i) = 0.50*std::log(rReturnMappingVariables.TrialEigenValues(i));
}







} // Namespace Kratos

