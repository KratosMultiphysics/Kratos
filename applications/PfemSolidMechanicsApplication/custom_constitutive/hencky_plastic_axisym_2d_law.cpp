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
#include "custom_constitutive/hencky_plastic_axisym_2d_law.hpp"

#include "solid_mechanics_application.h"

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

HenckyElasticPlasticAxisym2DLaw::HenckyElasticPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
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



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HenckyElasticPlasticAxisym2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
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

void HenckyElasticPlasticAxisym2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e= 0.5*(1-invbT*invb)
    Matrix InverseLeftCauchyGreen ( rLeftCauchyGreen.size1() , rLeftCauchyGreen.size2() );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 2, 2 ) );
    rStrainVector[3] = -InverseLeftCauchyGreen( 0, 1 );


}


//*********************** COMPUTE THE FIRST TERM OF THE CONSTITUTIVE MATRIX **********
//************************************************************************************

void HenckyElasticPlasticAxisym2DLaw::CalculateEigenValuesConstitutiveMatrix(const Matrix & rPrincipalTangent, const Matrix & rEigenVectors, Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();



    for(unsigned int i=0; i<4; i++)
    {
        for(unsigned int j=0; j<4; j++)
        {  
            rConstitutiveMatrix( i, j ) = EigenValuesConstitutiveComponent(rConstitutiveMatrix( i, j ), 
					  rPrincipalTangent, rEigenVectors, 
                                          this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
        }

    }

}

//*********************** COMPUTE THE SECOND TERM OF THE CONSTITUTIVE MATRIX**********
//************************************************************************************
void HenckyElasticPlasticAxisym2DLaw::CalculateEigenVectorsConstitutiveMatrix(const Matrix& rEigenVectors, const Vector& rEigenValues, const Matrix& rStressMatrix, Matrix& rConstitutiveMatrix)
{
    Vector PrincipalStress;
    PrincipalStress = this->GetStressVectorFromMatrix(rStressMatrix, PrincipalStress);
    
    rConstitutiveMatrix.clear();


    for(unsigned int i=0; i<4; i++)
    {
        for(unsigned int j=0; j<4; j++)
        {
            rConstitutiveMatrix( i, j ) = EigenVectorsConstitutiveComponent(rConstitutiveMatrix( i, j ), 
					  rEigenVectors, rEigenValues, PrincipalStress,
                                          this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
        }

    }

}











} // Namespace Kratos

