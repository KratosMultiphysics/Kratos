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
#include "custom_constitutive/hyperelastic_plastic_U_P_3D_law.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticUP3DLaw::HyperElasticPlasticUP3DLaw()
    : HyperElasticPlastic3DLaw()
{

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticUP3DLaw::HyperElasticPlasticUP3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
: HyperElasticPlastic3DLaw(pFlowRule,pYieldCriterion,pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticPlasticUP3DLaw::HyperElasticPlasticUP3DLaw(const HyperElasticPlasticUP3DLaw& rOther)
    : HyperElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticPlasticUP3DLaw::Clone() const
{
    return Kratos::make_shared<HyperElasticPlasticUP3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticPlasticUP3DLaw::~HyperElasticPlasticUP3DLaw()
{
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//******************************* COMPUTE DOMAIN PRESSURE  ***************************
//************************************************************************************


double &  HyperElasticPlasticUP3DLaw::CalculateVolumetricPressure (const MaterialResponseVariables & rElasticVariables,
							double & rPressure)
{

    const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  DomainGeometry.size();

    rPressure = 0;
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(PRESSURE);
    }

    return rPressure;

}

//************************* COMPUTE DOMAIN PRESSURE FACTORS***************************
//************************************************************************************

Vector&  HyperElasticPlasticUP3DLaw::CalculateVolumetricPressureFactors (const MaterialResponseVariables & rElasticVariables,
							      Vector & rFactors)

{
    double Pressure = 0;
    Pressure = this->CalculateVolumetricPressure( rElasticVariables, Pressure );

    if(rFactors.size()!=3) rFactors.resize(3);

    rFactors[0] =  1.0;
    rFactors[1] =  2.0;
    rFactors[2] =  Pressure*rElasticVariables.DeterminantF;

    return rFactors;
}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticPlasticUP3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
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
