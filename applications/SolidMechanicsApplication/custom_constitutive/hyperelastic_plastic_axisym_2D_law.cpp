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
#include "custom_constitutive/hyperelastic_plastic_axisym_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticAxisym2DLaw::HyperElasticPlasticAxisym2DLaw()
    : HyperElasticPlastic3DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticAxisym2DLaw::HyperElasticPlasticAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HyperElasticPlastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlasticAxisym2DLaw::HyperElasticPlasticAxisym2DLaw(const HyperElasticPlasticAxisym2DLaw& rOther)
    : HyperElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlasticAxisym2DLaw::Clone() const
{
    return Kratos::make_shared<HyperElasticPlasticAxisym2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticAxisym2DLaw::~HyperElasticPlasticAxisym2DLaw()
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

void HyperElasticPlasticAxisym2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
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

void HyperElasticPlasticAxisym2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
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

void HyperElasticPlasticAxisym2DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
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


void HyperElasticPlasticAxisym2DLaw::CalculateVolumetricConstitutiveMatrix ( const MaterialResponseVariables & rElasticVariables,
								  Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Vector Factors(3);
    noalias(Factors) = ZeroVector(3);
    Factors =this-> CalculateVolumetricPressureFactors( rElasticVariables, Factors );

    for(unsigned int i=0; i<3; i++)
    {
        for(unsigned int j=0; j<3; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, Factors,
                                          this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
        }

    }


}


//***********************COMPUTE PLASTIC CONSTITUTIVE MATRIX**************************
//************************************************************************************

void HyperElasticPlasticAxisym2DLaw::CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
									      FlowRule::RadialReturnVariables & rReturnMappingVariables,
									      Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Matrix IsoStressMatrix = rReturnMappingVariables.TrialIsoStressMatrix;

    //std::cout<< " TrialStressMatrix 2D "<<IsoStressMatrix<<std::endl;

    FlowRule::PlasticFactors ScalingFactors;
    mpFlowRule->CalculateScalingFactors( rReturnMappingVariables, ScalingFactors );

    for(unsigned int i=0; i<4; i++)
    {
        for(unsigned int j=0; j<4; j++)
        {
		rConstitutiveMatrix( i, j ) = PlasticConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, IsoStressMatrix, ScalingFactors,
                                          this->msIndexVoigt2D4C[i][0], this->msIndexVoigt2D4C[i][1], this->msIndexVoigt2D4C[j][0], this->msIndexVoigt2D4C[j][1]);
        }

    }


}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticPlasticAxisym2DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
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
