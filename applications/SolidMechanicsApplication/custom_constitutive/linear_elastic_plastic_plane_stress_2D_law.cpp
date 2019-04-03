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
#include "custom_constitutive/linear_elastic_plastic_plane_stress_2D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlasticPlaneStress2DLaw::LinearElasticPlasticPlaneStress2DLaw()
    : LinearElasticPlasticPlaneStrain2DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlasticPlaneStress2DLaw::LinearElasticPlasticPlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
 : LinearElasticPlasticPlaneStrain2DLaw(pFlowRule,pYieldCriterion,pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearElasticPlasticPlaneStress2DLaw::LinearElasticPlasticPlaneStress2DLaw(const LinearElasticPlasticPlaneStress2DLaw& rOther)
    : LinearElasticPlasticPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearElasticPlasticPlaneStress2DLaw::Clone() const
{
    return Kratos::make_shared<LinearElasticPlasticPlaneStress2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearElasticPlasticPlaneStress2DLaw::~LinearElasticPlasticPlaneStress2DLaw()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//********************* COMPUTE LINEAR ELASTIC CONSTITUTIVE MATRIX *******************
//************************************************************************************

void LinearElasticPlasticPlaneStress2DLaw::CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,
        const double& YoungModulus,
        const double& PoissonCoefficient )
{
    rLinearElasticMatrix.clear();

    // Plane stress constitutive matrix
    rLinearElasticMatrix ( 0 , 0 ) = (YoungModulus)/(1.0-PoissonCoefficient*PoissonCoefficient);
    rLinearElasticMatrix ( 1 , 1 ) = rLinearElasticMatrix ( 0 , 0 );

    rLinearElasticMatrix ( 2 , 2 ) = rLinearElasticMatrix ( 0 , 0 )*(1.0-PoissonCoefficient)*0.5;

    rLinearElasticMatrix ( 0 , 1 ) = rLinearElasticMatrix ( 0 , 0 )*PoissonCoefficient;
    rLinearElasticMatrix ( 1 , 0 ) = rLinearElasticMatrix ( 0 , 1 );
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearElasticPlasticPlaneStress2DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( PLANE_STRESS_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}

} // Namespace Kratos
